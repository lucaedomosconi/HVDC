/*!
 * @file
 * @author  Alessandro Lombardi
 * @version 0.1
 *
 * @section LICENSE
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * https://www.gnu.org/copyleft/gpl.html
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
#include <quad_operators.h>

//#include<test_1.h>
//#include<test_2.h>
//#include<test_3.h>
#include<test_4bis.h>
//#include<test_5.h>

int
main (int argc, char **argv)
{

  using q1_vec  = q1_vec<distributed_vector>;
  
  /*
  Manegement of solutions ordering: ord0-> phi                 ord1->rho
  Equation ordering: ord0->diffusion-reaction equation         ord1->continuity equation 
  */
  ordering
    ord0 = [] (tmesh::idx_t gt) -> size_t { return dof_ordering<2, 0> (gt); },
    ord1 = [] (tmesh::idx_t gt) -> size_t { return dof_ordering<2, 1> (gt); };
    
  // Initialize MPI
  MPI_Init (&argc, &argv);
  int rank, size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  // Generate the mesh in 3d
  tmesh tmsh;
  tmsh.read_connectivity (simple_conn_p, simple_conn_num_vertices,
                          simple_conn_t, simple_conn_num_trees);

  //Uniform refinement
  int recursive = 1;
  tmsh.set_refine_marker (uniform_refinement);
  tmsh.refine (recursive);

  //In test 1 we only have uniform refinement, in all other cases we perform additional refinement
  if (extra_refinement) 
  {
    tmsh.set_refine_marker(refinement);
    tmsh.refine (recursive);

    tmsh.set_coarsen_marker(coarsening);
    tmsh.coarsen(recursive);
  }

  tmesh::idx_t gn_nodes = tmsh.num_global_nodes ();
  tmesh::idx_t ln_nodes = tmsh.num_owned_nodes ();
  tmesh::idx_t ln_elements = tmsh.num_local_quadrants (); 

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
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep ();
       ++quadrant)
    {
      double xx{quadrant->centroid(0)}, yy{quadrant->centroid(1)};  
      
      epsilon[quadrant->get_forest_quad_idx ()] = epsilon_fun(xx,yy);
      sigma[quadrant->get_forest_quad_idx ()] = sigma_fun(xx,yy);
      
      delta0[quadrant->get_forest_quad_idx ()] = -1.0;
      delta1[quadrant->get_forest_quad_idx ()] = 1.0;
      f0[quadrant->get_forest_quad_idx ()] = 0.0;
      f1[quadrant->get_forest_quad_idx ()] = 1.0;

      for (int ii = 0; ii < 4; ++ii)
        {
          if (! quadrant->is_hanging (ii))
          {
            zero[quadrant->gt (ii)] = 0.;
            zeta0[quadrant->gt (ii)] = 1.0;
            zeta1[quadrant->gt (ii)] = 1.0;
            g0[quadrant->gt (ii)] = 0.;

            sold[ord0(quadrant->gt (ii))] = 0.0;
            sold[ord1(quadrant->gt (ii))] = 0.0;
            sol[ord0(quadrant->gt (ii))] = 0.0;
            sol[ord1(quadrant->gt (ii))] = 0.0;
          }
          else
            for (int jj = 0; jj < 2; ++jj)
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
  tmsh.octbin_export_quadrant ("epsilon_file", epsilon);
  tmsh.octbin_export_quadrant ("sigma_file", sigma);

  bim2a_solution_with_ghosts (tmsh, sold, replace_op, ord0, false);
  bim2a_solution_with_ghosts (tmsh, sold, replace_op, ord1);

  zero.assemble (replace_op);
  zeta0.assemble (replace_op);
  zeta1.assemble (replace_op);
  g0.assemble (replace_op);                      

  // Save inital conditions
  sprintf(filename, "model_0_u_0000");
  tmsh.octbin_export (filename, sold, ord0);
  sprintf(filename, "model_0_v_0000");
  tmsh.octbin_export (filename, sold, ord1); 

  int count = 0;

  //Choosing the indices for the nodes corresponding to the vales of rho of interest 
  if (rank == 0) {
    rho_out[0]=std::vector<double>(N_rhos+1,0.0);
    rho_idx = find_idx(tmsh,points,tols,N_rhos);
  }
 
  // Time cycle
  for( double time = DELTAT; time <= T; time += DELTAT)
  {
    count++;

    // Define boundary conditions
    dirichlet_bcs bcs0, bcs1;
    bcs0.push_back (std::make_tuple (0, 2, [](double x, double y){return 0.0;})); //bottom
    bcs0.push_back (std::make_tuple (0, 3, [time](double x, double y){return 1.5e4 * (1 - exp(-time/tau));})); //top

    // Print curent time
    if(rank==0)
      std::cout<<"TIME= "<<time<<std::endl;

    // Reset containers
    A.reset ();
    sol.get_owned_data ().assign (sol.get_owned_data ().size (), 0.0);
    sol.assemble (replace_op);

    // Initialize non constant (in time) parameters
    for (auto quadrant = tmsh.begin_quadrant_sweep ();
    quadrant != tmsh.end_quadrant_sweep ();
    ++quadrant)
    {
      for (int ii = 0; ii < 4; ++ii)
        if (! quadrant->is_hanging (ii))
          g1[quadrant->gt (ii)] = sold[ord1(quadrant->gt (ii))];
          
        else
          for (int jj = 0; jj < 2; ++jj)
            g1[quadrant->gparent (jj, ii)] += 0.;
    }

    g1.assemble(replace_op);

    // advection_diffusion
    bim2a_advection_diffusion (tmsh, epsilon, zero, A, true, ord0, ord0);
    bim2a_advection_diffusion (tmsh, sigma, zero, A, true, ord1, ord0);

    // reaction
    bim2a_reaction (tmsh, delta0, zeta0, A, ord0, ord1);
    bim2a_reaction (tmsh, delta1, zeta1, A, ord1, ord1);

    //rhs
    bim2a_rhs (tmsh, f0, g0, sol, ord0);
    bim2a_rhs (tmsh, f1, g1, sol, ord1);

    //boundary conditions
    bim2a_dirichlet_bc (tmsh, bcs0, A, sol, ord1, ord0, false);

    // Communicate matrix and RHS
    A.assemble ();
    sol.assemble ();

    // Solver analysis
    lin_solver->set_lhs_distributed ();
    A.aij (xa, ir, jc, lin_solver->get_index_base ());
    lin_solver->set_distributed_lhs_structure (A.rows (), ir, jc);
    std::cout << "lin_solver->analyze () return value = "<< lin_solver->analyze () << std::endl;

    // Matrix update
    A.aij_update (xa, ir, jc, lin_solver->get_index_base ());
    lin_solver->set_distributed_lhs_data (xa);

    // Factorization
    std::cout << "lin_solver->factorize () = " << lin_solver->factorize () << std::endl;

    // Set RHS data
    lin_solver->set_rhs_distributed (sol);

    // Solution
    std::cout << "lin_solver->solve () = " << lin_solver->solve () << std::endl;

    // Copy solution
    q1_vec result = lin_solver->get_distributed_solution ();
    for (int idx = sold.get_range_start (); idx < sold.get_range_end (); ++idx)
      sold (idx) = result (idx);
    sold.assemble (replace_op);
    
    // Save solution
    if (save_sol == true)
    {
      sprintf(filename, "model_0_u_%4.4d",count);
      tmsh.octbin_export (filename, sold, ord0);
      sprintf(filename, "model_0_v_%4.4d",count);
      tmsh.octbin_export (filename, sold, ord1);      
    }

    // Save rho values
    if (rank == 0)
    {
      std::vector<double> temp(N_rhos+1);
      temp[0] = time;
      for (size_t i=1; i < N_rhos+1; i++)
        temp[i] = sold[ord1(rho_idx[i-1])];

      rho_out[count] = temp;
    }
  }
  
  // Print file with rho values
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

  // Clean linear solver
  lin_solver->cleanup ();
  
  MPI_Finalize (); 

  return 0;
}
