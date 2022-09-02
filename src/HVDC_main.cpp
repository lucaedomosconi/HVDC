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
#include <tmesh_3d.h>
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

  // Manegement of solutions ordering, ord0-> phi    ord1->rho
  // Equation ordering: ord0->other equation         ord1->continuity equation 
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
  mumps *lin_solver2 = new mumps ();

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

  //bool find_idx = true;
  //int rho_idx{1};
  
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
          /*if (find_idx && std::fabs(quadrant->p(2, ii) - 0.000984375) < 1e-9 && rank == 0)
          {
            find_idx = false;
            rho_idx = quadrant -> t(ii);
          }*/
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
  //bim3a_solution_with_ghosts (tmsh, zero, replace_op);
  //bim3a_solution_with_ghosts (tmsh, zeta, replace_op);
  //bim3a_solution_with_ghosts (tmsh, g, replace_op);

  bim3a_solution_with_ghosts (tmsh, sold, replace_op, ord0, false);
  bim3a_solution_with_ghosts (tmsh, sold, replace_op, ord1);

  zero.assemble (replace_op);
  zeta0.assemble (replace_op);
  zeta1.assemble (replace_op);
  g0.assemble (replace_op);         
  //sold.assemble(replace_op);              
  TOC ("compute coefficient");

  // Save inital conditions
  sprintf(filename, "model_0_u_0000");
  tmsh.octbin_export (filename, sold, ord0);
  sprintf(filename, "model_0_v_0000");
  tmsh.octbin_export (filename, sold, ord1); 
  
  
  // advection_diffusion
  TIC ();
  bim3a_advection_diffusion (tmsh, epsilon, zero, A, true, ord0, ord0);
  bim3a_advection_diffusion (tmsh, sigma, zero, A, true, ord1, ord0);

  // reaction
  bim3a_reaction (tmsh, delta0, zeta0, A, ord0, ord1);
  bim3a_reaction (tmsh, delta1, zeta1, A, ord1, ord1);
  TOC ("assemble LHS");
  
  TIC();
  //bim3a_rhs (tmsh, f0, g0, sol, ord0);
  //bim3a_rhs (tmsh, f1, g0, sol, ord1); 
  TOC ("assemble RHS");

  //boundary conditions
  TIC ();
  //bim3a_dirichlet_bc (tmsh, bcs0, A, sol, ord0, false);
  //bim3a_dirichlet_bc (tmsh, bcs1, A, sol, ord1, false);
  TOC ("apply BCS");
  
  // Communicate matrix and RHS
  TIC ();
  A.assemble ();
  sol.assemble ();
  TOC ("communicate A");
  
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


  int count = 0;

  //std::cout << sold[ord0(5)] << std::endl;
  std::vector<std::pair<double,double>> rho_out;
  if (rank == 0)
    rho_out.push_back(std::make_pair(0.0,0.0));

  // Time cycle
  for( double time = DELTAT; time <= T; time += DELTAT)
  {
    count++;

    dirichlet_bcs3 bcs0, bcs1;
    
    bcs0.push_back (std::make_tuple (0, 4, [](double x, double y, double z){return 0.0;})); //bottom
    bcs0.push_back (std::make_tuple (0, 5, [time](double x, double y, double z){return 1.5e4 * (1 - exp(-time/tau));})); //top
    //bcs1.push_back (std::make_tuple (0, 4, [](double x, double y, double z){return 0.0;})); //bottom
    //bcs1.push_back (std::make_tuple (0, 5, [](double x, double y, double z){return 0.0;})); //top

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
    //bim3a_dirichlet_bc (tmsh, bcs1, A, sol, ord1, false);
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

    //Rho output
    /*func3_quad fun_one = [] (tmesh_3d::quadrant_iterator, tmesh_3d::idx_t) {return 1;};

    q1_vec Area(ln_nodes*2), charge(ln_nodes*2);
    charge.get_owned_data ().assign (charge.get_owned_data ().size (), 0.0);
    Area.get_owned_data ().assign (Area.get_owned_data ().size (), 0.0);
    
    bim3a_boundary_mass(tmsh, 0, 4, Area, fun_one, ord1);
    
    for (int idx = sold.get_range_start (); idx < sold.get_range_end (); ++idx)
    {
      charge(idx) = sold(idx) * Area(idx) * 1e9;
    }
    charge.assemble();*/
  
    if (rank == 0)
    {
      //double Q{0.0};
      //for (auto &i: charge.get_owned_data())
       // Q += i;
    
      //rho_out.push_back(std::make_pair(time, Q));
      rho_out.push_back(std::make_pair(time, sold[ord1(0)] ));
      //std::cout << "rho= " << sold[ord1(0)] << std::endl;
    }
  }
  
  if (rank == 0)
  {
    std::ofstream outFile("Arho.txt");
    for (const auto &e : rho_out)
      outFile << e.first << "\t" << e.second << "\n";
  }
 /*
 if(rank == 0)
  {
    octave_io_mode m = gz_write_mode;  
    // save data to file
    A.aij (xa, ir, jc, lin_solver->get_index_base ());
    ColumnVector oct_xa (xa.size (), 0.0);
    Array<octave_idx_type> oct_ir (dim_vector (ir.size (), 1), 0);
    Array<octave_idx_type> oct_jc (dim_vector (jc.size (), 1), 0);

    std::copy_n (xa.begin (), xa.size (), oct_xa.fortran_vec ());
    std::copy_n (ir.begin (), ir.size (), oct_ir.fortran_vec ());
    std::copy_n (jc.begin (), jc.size (), oct_jc.fortran_vec ());

    octave_scalar_map the_map;
    the_map.assign ("xa", oct_xa);
    the_map.assign ("ir", oct_ir);
    the_map.assign ("jc", oct_jc);

    ColumnVector oct_c (Jc.size(), 0.0), oct_d (Jd.size(), 0.0), oct_s (sold.size(), 0.0), oct_l (lambda.size(), 0.0);
    std::copy_n (Jc.get_owned_data().begin (), Jc.get_owned_data().size (), oct_c.fortran_vec ());
    std::copy_n (Jd.get_owned_data().begin (), Jd.get_owned_data().size (), oct_d.fortran_vec ());
    std::copy_n (sold.get_owned_data().begin (), sold.get_owned_data().size (), oct_s.fortran_vec ());
    std::copy_n (lambda.get_owned_data().begin (), lambda.get_owned_data().size (), oct_l.fortran_vec ());

    //octave_scalar_map the_map;
    the_map.assign ("Jd", oct_d);
    the_map.assign ("Jc", oct_c);
    the_map.assign ("sold", oct_s);
    the_map.assign ("lambda", oct_l);

    assert (octave_io_open ("currents.octbin", m, &m) == 0);
    assert (octave_save ("the_map", octave_value (the_map)) == 0);
    assert (octave_io_close () == 0);
  }
*/  
  // Close MPI and print report
  MPI_Barrier (MPI_COMM_WORLD);
  if (rank == 0) { print_timing_report (); }

  //if (rank == 0)
  //  std::cout << "Tau = " << (epsilon_0*epsilon_r)/sigma_ << std::endl;

  // Clean linear solver
  lin_solver->cleanup ();
  
  MPI_Finalize (); 

  return 0;
}