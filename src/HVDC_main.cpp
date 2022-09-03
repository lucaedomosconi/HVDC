/*
  Copyright (C) 2019 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <octave_file_io.h>

#include <bim_distributed_vector.h>
#include <bim_sparse_distributed.h>
#include <bim_timing.h>
#include <mumps_class.h>
#include <quad_operators_3d.h>

//#include<test_1.h>
#include<test_2.h>
//#include<test_3.h>
//#include<test_4.h>
//#include<test_5.h>


#define TIC() MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
#define TOC(S) MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { toc (S); }

int
main (int argc, char **argv)
{

  using q1_vec  = q1_vec<distributed_vector>;
  /*!
  Manegement of solutions ordering: ord0-> phi                 ord1->rho
  Equation ordering: ord0->diffusion-reaction equation         ord1->continuity equation 
  */
  ordering
    ord0 = [] (tmesh_3d::idx_t gt) -> size_t { return dof_ordering<2, 0> (gt); },
    ord1 = [] (tmesh_3d::idx_t gt) -> size_t { return dof_ordering<2, 1> (gt); };
    
  // Initialize MPI
  MPI_Init (&argc, &argv);
  int rank, size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  // Generate the mesh in 3d
  tmesh_3d tmsh;
  tmsh.read_connectivity (simple_conn_p, simple_conn_num_vertices,
                          simple_conn_t, simple_conn_num_trees);


  int recursive = 1;
  tmsh.set_refine_marker (uniform_refinement);
  tmsh.refine (recursive);

  /*!
  Only in test 1 we have uniform refinement, in all other cases we perform additional refinement
  */
  if (extra_refinement) 
  {
    tmsh.set_refine_marker(refinement);
    tmsh.refine (recursive);

    tmsh.set_coarsen_marker(coarsening);
    tmsh.coarsen(recursive);
  }

  tmesh_3d::idx_t gn_nodes = tmsh.num_global_nodes ();
  tmesh_3d::idx_t ln_nodes = tmsh.num_owned_nodes ();
  tmesh_3d::idx_t ln_elements = tmsh.num_local_quadrants (); 

  // Allocate linear solver
  mumps *lin_solver = new mumps ();

 // Allocate initial data container
  q1_vec sold (ln_nodes * 2);
  sold.get_owned_data ().assign (sold.get_owned_data ().size (), 0.0);

  q1_vec sol (ln_nodes * 2);
  sol.get_owned_data ().assign (sol.get_owned_data ().size (), 0.0);

  std::vector<double> xa;
  std::vector<int> ir, jc;
  
  // Declare system matrix
  distributed_sparse_matrix A;
  A.set_ranges (ln_nodes * 2);
  
  // Buffer for export filename
  char filename[255]="";

  //Output rho vector
  size_t N_timesteps = (size_t) (ceil(T/DELTAT)+1);
  std::vector<std::vector<double>> rho_out(N_timesteps,std::vector<double>(N_rhos+1));

  // Compute coefficients

  // diffusion
  std::vector<double> epsilon (ln_elements, 0.);
  std::vector<double> sigma (ln_elements, 0.);
  q1_vec zero (ln_nodes);

  // reaction
  std::vector<double> delta0 (ln_elements, 0.);
  std::vector<double> delta1 (ln_elements, 0.);
  q1_vec zeta0 (ln_nodes);
  q1_vec zeta1 (ln_nodes);

  // rhs
  std::vector<double> f0 (ln_elements, 0.);
  std::vector<double> f1 (ln_elements, 0.);
  q1_vec g0 (ln_nodes);
  q1_vec g1 (ln_nodes);

  // Initialize constant (in time) parameters and initial data
  TIC ();
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep ();
       ++quadrant)
    {
      double xx{quadrant->centroid(0)}, yy{quadrant->centroid(1)}, zz{quadrant->centroid(2)};  
      
      epsilon[quadrant->get_forest_quad_idx ()] = epsilon_fun(xx,yy,zz);
      sigma[quadrant->get_forest_quad_idx ()] = sigma_fun(xx,yy,zz);
      
      delta0[quadrant->get_forest_quad_idx ()] = -1.0;
      delta1[quadrant->get_forest_quad_idx ()] = 1.0;
      f0[quadrant->get_forest_quad_idx ()] = 0.0;
      f1[quadrant->get_forest_quad_idx ()] = 1.0;

      for (int ii = 0; ii < 8; ++ii)
        {
          if (! quadrant->is_hanging (ii))
          {
            zero[quadrant->gt (ii)] = 0.;
            zeta0[quadrant->gt (ii)] = 1.0;
            zeta1[quadrant->gt (ii)] = 1.0;
            g0[quadrant->gt (ii)] = 0.;

            double zz=quadrant->p(2,ii);
            sold[ord0(quadrant->gt (ii))] = 0.0;
            sold[ord1(quadrant->gt (ii))] = 0.0;
            sol[ord0(quadrant->gt (ii))] = 0.0;
            sol[ord1(quadrant->gt (ii))] = 0.0;
          }
          else
            for (int jj = 0; jj < quadrant->num_parents (ii); ++jj)
              {
                zero[quadrant->gparent (jj, ii)] += 0.;
                zeta0[quadrant->gparent (jj, ii)] += 0.;
                zeta1[quadrant->gparent (jj, ii)] += 0.;
                g0[quadrant->gparent (jj, ii)] += 0.;

                sold[ord0(quadrant->gparent (jj, ii))] += 0.;
                sold[ord1(quadrant->gparent (jj, ii))] += 0.;                
              }
        }
    }
  
  bim3a_solution_with_ghosts (tmsh, sold, replace_op, ord0, false);
  bim3a_solution_with_ghosts (tmsh, sold, replace_op, ord1);

  zero.assemble (replace_op);
  zeta0.assemble (replace_op);
  zeta1.assemble (replace_op);
  g0.assemble (replace_op);                      
  TOC ("compute coefficient");

  // Save inital conditions
  sprintf(filename, "model_0_u_0000");
  tmsh.octbin_export (filename, sold, ord0);
  sprintf(filename, "model_0_v_0000");
  tmsh.octbin_export (filename, sold, ord1); 

  int count = 0;
 
  if (rank == 0) {
    rho_out[0]=std::vector<double>(N_rhos+1,0.0);
    rho_idx = find_idx(tmsh,points,tols,N_rhos);
  }
 
  // Time cycle
  for( double time = DELTAT; time <= T; time += DELTAT)
  {
    count++;

    dirichlet_bcs3 bcs0, bcs1;
    
    bcs0.push_back (std::make_tuple (0, 4, [](double x, double y, double z){return 0.0;})); //bottom
    bcs0.push_back (std::make_tuple (0, 5, [time](double x, double y, double z){return 1.5e4 * (1 - exp(-time/tau));})); //top

    // Print curent time
    if(rank==0)
      std::cout<<"TIME= "<<time<<std::endl;

    // Reset containers
    TIC();
    A.reset ();
    sol.get_owned_data ().assign (sol.get_owned_data ().size (), 0.0);
    sol.assemble (replace_op);
    TOC("Resetting")

    // Initialize non constant (in time) parameters
    TIC();
    for (auto quadrant = tmsh.begin_quadrant_sweep ();
    quadrant != tmsh.end_quadrant_sweep ();
    ++quadrant)
    {
      for (int ii = 0; ii < 8; ++ii)
        if (! quadrant->is_hanging (ii))
          g1[quadrant->gt (ii)] = sold[ord1(quadrant->gt (ii))];
          
        else
          for (int jj = 0; jj < quadrant->num_parents (ii); ++jj)
            g1[quadrant->gparent (jj, ii)] += 0.;
    }

    g1.assemble(replace_op);
    TOC("Update rhs");

    // advection_diffusion
    TIC ();
    bim3a_advection_diffusion (tmsh, epsilon, zero, A, true, ord0, ord0);
    bim3a_advection_diffusion (tmsh, sigma, zero, A, true, ord1, ord0);

    // reaction
    bim3a_reaction (tmsh, delta0, zeta0, A, ord0, ord1);
    bim3a_reaction (tmsh, delta1, zeta1, A, ord1, ord1);
    TOC ("assemble LHS");

    //rhs
    TIC ();
    bim3a_rhs (tmsh, f0, g0, sol, ord0);
    bim3a_rhs (tmsh, f1, g1, sol, ord1);
    TOC ("assemble RHS");

    //boundary conditions
    TIC ();
    bim3a_dirichlet_bc (tmsh, bcs0, A, sol, ord1, ord0, false);
    TOC ("apply BCS");

    // Communicate matrix and RHS
    TIC ();
    A.assemble ();
    sol.assemble ();
    TOC ("communicate A & b");

    // Solver analysis
    TIC ();
    lin_solver->set_lhs_distributed ();
    A.aij (xa, ir, jc, lin_solver->get_index_base ());
    lin_solver->set_distributed_lhs_structure (A.rows (), ir, jc);
    std::cout << "lin_solver->analyze () return value = "<< lin_solver->analyze () << std::endl;
    TOC ("solver analysis");

    // Matrix update
    TIC ();
    A.aij_update (xa, ir, jc, lin_solver->get_index_base ());
    lin_solver->set_distributed_lhs_data (xa);
    TOC ("set LHS data");

    // Factorization
    TIC ();
    std::cout << "lin_solver->factorize () = " << lin_solver->factorize () << std::endl;
    TOC ("solver factorize");

    // Set RHS data
    TIC ();
    lin_solver->set_rhs_distributed (sol);
    TOC ("set RHS data");

    // Solution
    TIC ();
    std::cout << "lin_solver->solve () = " << lin_solver->solve () << std::endl;
    TOC ("solver solve");

    // Copy solution
    TIC();
    q1_vec result = lin_solver->get_distributed_solution ();
    for (int idx = sold.get_range_start (); idx < sold.get_range_end (); ++idx)
      sold (idx) = result (idx);
    sold.assemble (replace_op);
    TOC("Obtaining solution");
    
    // Save solution
    if (save_sol == true)
    {
      TIC();
      sprintf(filename, "model_0_u_%4.4d",count);
      tmsh.octbin_export (filename, sold, ord0);
      sprintf(filename, "model_0_v_%4.4d",count);
      tmsh.octbin_export (filename, sold, ord1);      
      TOC("Exporting solution");
    }
  
    if (rank == 0)
    {
      std::vector<double> temp(N_rhos+1);
      temp[0] = time;
      for (size_t i=1; i < N_rhos+1; i++)
        temp[i] = sold[ord1(rho_idx[i-1])];

      rho_out[count] = temp;
    }
  }
  
  if (rank == 0)
  {
    std::ofstream outFile("Arho.txt");
    outFile << N_rhos+1 << "\t" << N_timesteps << "\n";

    for (const auto &e : rho_out) {
      for (size_t i=0; i < N_rhos+1; i++)
        outFile << e[i] << "\t";
      outFile << "\n";
    }
  }
  
  // Close MPI and print report
  MPI_Barrier (MPI_COMM_WORLD);
  if (rank == 0) { print_timing_report (); }

  // Clean linear solver
  lin_solver->cleanup ();
  
  MPI_Finalize (); 

  return 0;
}