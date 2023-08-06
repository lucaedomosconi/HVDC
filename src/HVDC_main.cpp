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
#include <iomanip>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <octave_file_io.h>
#include <array>
#include <functional>

#include <bim_distributed_vector.h>
#include <bim_sparse_distributed.h>
#include <bim_timing.h>
#include <mumps_class.h>
#include <quad_operators.h>

#include <nlohmann/json.hpp>
#include<test_1_2d.h>
//#include<test_2.h>
//#include<test_3.h>
//#include<test_4bis.h>
//#include<test_5.h>

constexpr size_t N_eqs= 5;
/*
template<size_t... args>
std::array<ordering,N_eqs> makeorder(){
  return std::array<ordering,N_eqs> {dof_ordering<N_eqs,args>...};
}
*//*
template<class>
std::array<ordering,N_eqs> makeorder_impl();

template<size_t... args>
std::array<ordering,N_eqs> makeorder_impl<std::integer_sequence<size_t,args...>>(){
  return std::array<ordering,N_eqs> {dof_ordering<N_eqs,args>...};
};

template<size_t N>
auto makeorder(){
  return makeorder_impl<std::make_integer_sequence<size_t,N>>();
}
*/

template<class>
struct make_order_struct{};

template<size_t... args>
struct make_order_struct<std::integer_sequence<size_t,args...>> {
  static std::array<ordering,N_eqs> fun(){
    return std::array<ordering,N_eqs> {dof_ordering<N_eqs,args>...};
  }
};

template<size_t N>
auto makeorder(){
  return make_order_struct<std::make_integer_sequence<size_t,N>>::fun();
}

using q1_vec_ = q1_vec<distributed_vector>;

void time_step(const int rank, const double time, const double DELTAT,
                const std::array<ordering,N_eqs> &ord,
                tmesh &tmsh, mumps *lin_solver, distributed_sparse_matrix &A,
                std::vector<double> &xa, std::vector<int> &ir, std::vector<int> &jc,
                std::vector<double> &epsilon, std::vector<double> &sigma,
                std::vector<double> &zero_std_vect, q1_vec_ &zero_q1,
                std::vector<double> &delta1, std::vector<double> &delta0,
                std::vector<double> &reaction_term_p1,
                std::vector<double> &reaction_term_p2,
                std::vector<double> &reaction_term_p3,
                std::vector<double> &diffusion_term_p1,
                std::vector<double> &diffusion_term_p2,
                std::vector<double> &diffusion_term_p3,
                q1_vec_ &zeta0, q1_vec_ &zeta1,
                std::vector<double> &f1, std::vector<double> &f0,
                q1_vec_ &g1, q1_vec_ &g0, q1_vec_ &gp1, q1_vec_ &gp2, q1_vec_ &gp3,
                q1_vec_ &sold, q1_vec_ &sol) {

    // Define boundary conditions
    dirichlet_bcs bcs0, bcs1;
    bcs0.push_back (std::make_tuple (0, 0, [](double x, double y){return 0.0;})); //bottom
    bcs0.push_back (std::make_tuple (0, 1, [time,DELTAT](double x, double y){return 1.5e4 * (1 - exp(-(time+DELTAT)/tau));})); //top

    // Print curent time
    if(rank==0)
      std::cout << "TIME= "<< time + DELTAT <<std::endl;

    // Reset containers
    A.reset ();
    sol.get_owned_data ().assign (sol.get_owned_data ().size (), 0.0);
    sol.assemble (replace_op);

    // Initialize non constant (in time) parameters
    for (auto quadrant = tmsh.begin_quadrant_sweep ();
    quadrant != tmsh.end_quadrant_sweep ();
    ++quadrant)
    {
      double xx{quadrant->centroid(0)}, yy{quadrant->centroid(1)};  

      reaction_term_p1[quadrant->get_forest_quad_idx ()] =
        1 + DELTAT / tau_p1;
      reaction_term_p2[quadrant->get_forest_quad_idx ()] =
        1 + DELTAT / tau_p2;
      reaction_term_p3[quadrant->get_forest_quad_idx ()] =
        1 + DELTAT / tau_p3;
      diffusion_term_p1[quadrant->get_forest_quad_idx ()] =
        - DELTAT / tau_p1 * epsilon_0 * csi_1_fun(xx,yy);
      diffusion_term_p2[quadrant->get_forest_quad_idx ()] =
        - DELTAT / tau_p2 * epsilon_0 * csi_2_fun(xx,yy);
      diffusion_term_p3[quadrant->get_forest_quad_idx ()] =
        - DELTAT / tau_p3 * epsilon_0 * csi_3_fun(xx,yy);
      sigma[quadrant->get_forest_quad_idx ()] = sigma_fun(xx,yy,DELTAT);
      for (int ii = 0; ii < 4; ++ii)
        if (! quadrant->is_hanging (ii)){
          g0[quadrant->gt (ii)] = sold[ord[0](quadrant->gt (ii))];

          gp1[quadrant->gt (ii)] = sold[ord[2](quadrant->gt (ii))];
          gp2[quadrant->gt (ii)] = sold[ord[3](quadrant->gt (ii))];
          gp3[quadrant->gt (ii)] = sold[ord[4](quadrant->gt (ii))];
        }
        else
          for (int jj = 0; jj < 2; ++jj){
            g0[quadrant->gparent (jj, ii)] += 0.;
            gp1[quadrant->gparent (jj, ii)] += 0.;
            gp2[quadrant->gparent (jj, ii)] += 0.;
            gp3[quadrant->gparent (jj, ii)] += 0.;
          }
    }

    g0.assemble(replace_op);
    gp1.assemble(replace_op);
    gp2.assemble(replace_op);
    gp3.assemble(replace_op);
    // advection_diffusion

    bim2a_advection_diffusion (tmsh, sigma, zero_q1, A, true, ord[0], ord[1]);
    bim2a_advection_diffusion (tmsh, epsilon, zero_q1, A, true, ord[1], ord[1]);
    bim2a_advection_diffusion (tmsh, diffusion_term_p1, zero_q1, A, true, ord[2], ord[1]);
    bim2a_advection_diffusion (tmsh, diffusion_term_p2, zero_q1, A, true, ord[3], ord[1]);
    bim2a_advection_diffusion (tmsh, diffusion_term_p3, zero_q1, A, true, ord[4], ord[1]);
    
    // reaction
    bim2a_reaction (tmsh, delta0, zeta0, A, ord[0], ord[0]);
    bim2a_reaction (tmsh, delta1, zeta1, A, ord[1], ord[0]);
    bim2a_reaction (tmsh, delta0, zeta0, A, ord[1], ord[2]);
    bim2a_reaction (tmsh, delta0, zeta0, A, ord[1], ord[3]);
    bim2a_reaction (tmsh, delta0, zeta0, A, ord[1], ord[4]);
    bim2a_reaction (tmsh, reaction_term_p1, zeta1, A, ord[2], ord[2]);
    bim2a_reaction (tmsh, reaction_term_p1, zeta1, A, ord[3], ord[3]);
    bim2a_reaction (tmsh, reaction_term_p1, zeta1, A, ord[4], ord[4]);

    //rhs
    bim2a_rhs (tmsh, f0, g0, sol, ord[0]);
    bim2a_rhs (tmsh, f1, g1, sol, ord[1]);
    bim2a_rhs (tmsh, f0, gp1, sol, ord[2]);
    bim2a_rhs (tmsh, f0, gp1, sol, ord[3]);
    bim2a_rhs (tmsh, f0, gp1, sol, ord[4]);

    //boundary conditions
    bim2a_dirichlet_bc (tmsh, bcs0, A, sol, ord[0], ord[1], false);

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
    q1_vec<distributed_vector> result = lin_solver->get_distributed_solution ();
    for (int idx = sold.get_range_start (); idx < sold.get_range_end (); ++idx)
      sold (idx) = result (idx);
    sold.assemble (replace_op);

    
    /*
    // Save rho values
    if (rank == 0)
    {
      std::vector<double> temp(N_rhos+1);
      temp[0] = time;
      for (size_t i=1; i < N_rhos+1; i++)
        temp[i] = sold[ord[0](rho_idx[i-1])];

      rho_out[count] = temp;
      
    }
    */
}






int
main (int argc, char **argv)
{

  using q1_vec = q1_vec<distributed_vector>;
  using json = nlohmann::json;

  std::ifstream data_file("data.json");
  json data = json::parse(data_file);

  T = data["T"];
  tau = data["tau"];
  tau_p1 = data["tau_p1"];
  tau_p2 = data["tau_p2"];
  tau_p3 = data["tau_p3"];
  save_sol = data["save_sol"];
  
  epsilon_0 = data["epsilon_0"];
  epsilon_inf = data["epsilon_inf"];           // permittivity at infinite frequency
  csi1 = data["csi1"];
  csi2 = data["csi2"];
  csi3 = data["csi3"];
  sigma_ = data["sigma_"];                     // conducivity coeff

  double DT = data["DT"];
  double toll = data["toll_of_adaptive_time_step"];
  bool save_error = data["save_error"];

  /*
  Manegement of solutions ordering: ord[0]-> phi                 ord[1]->rho
  Equation ordering: ord[0]->diffusion-reaction equation         ord[1]->continuity equation 
  */
  
//  const std::array<ordering,N_eqs> ord{dof_ordering<N_eqs,0>,
//                                      dof_ordering<N_eqs,1>};
  const std::array<ordering,N_eqs> ord(makeorder<N_eqs>());


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
  q1_vec sold (ln_nodes * N_eqs);
  sold.get_owned_data ().assign (sold.get_owned_data ().size (), 0.0);

  q1_vec sol (ln_nodes * N_eqs);
  sol.get_owned_data ().assign (sol.get_owned_data ().size (), 0.0);

  std::vector<double> xa;
  std::vector<int> ir, jc;
  
  // Declare system matrix
  distributed_sparse_matrix A;
  A.set_ranges (ln_nodes * N_eqs);
  
  // Buffer for export filename
  char filename[255]="";

  //Output rho vector
  //size_t N_timesteps = (size_t) (ceil(T/DELTAT)+1);
  //std::vector<std::vector<double>> rho_out(N_timesteps,std::vector<double>(N_rhos+1));

  // Compute coefficients

  // diffusion
  std::vector<double> epsilon (ln_elements, 0.);
  std::vector<double> sigma (ln_elements, 0.);
  std::vector<double> zero_std_vect(ln_elements, 0.);
  std::vector<double> diffusion_term_p1 (ln_elements,0.);
  std::vector<double> diffusion_term_p2 (ln_elements,0.);
  std::vector<double> diffusion_term_p3 (ln_elements,0.);
  q1_vec zero_q1 (ln_nodes);

  // reaction
  std::vector<double> delta1 (ln_elements, 0.);
  std::vector<double> delta0 (ln_elements, 0.);
  std::vector<double> reaction_term_p1 (ln_elements,0.);
  std::vector<double> reaction_term_p2 (ln_elements,0.);
  std::vector<double> reaction_term_p3 (ln_elements,0.);
  q1_vec zeta0 (ln_nodes);
  q1_vec zeta1 (ln_nodes);

  // rhs
  std::vector<double> f1 (ln_elements, 0.);
  std::vector<double> f0 (ln_elements, 0.);
  q1_vec g1 (ln_nodes);
  q1_vec g0 (ln_nodes);
  q1_vec gp1 (ln_nodes);
  q1_vec gp2 (ln_nodes);
  q1_vec gp3 (ln_nodes);

  // Initialize constant (in time) parameters and initial data
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep ();
       ++quadrant)
    {
      double xx{quadrant->centroid(0)}, yy{quadrant->centroid(1)};  
      
      epsilon[quadrant->get_forest_quad_idx ()] = epsilon_fun(xx,yy);
      
      delta1[quadrant->get_forest_quad_idx ()] = -1.0;
      delta0[quadrant->get_forest_quad_idx ()] = 1.0;
      f1[quadrant->get_forest_quad_idx ()] = 0.0;
      f0[quadrant->get_forest_quad_idx ()] = 1.0;

      for (int ii = 0; ii < 4; ++ii)
        {
          if (! quadrant->is_hanging (ii))
          {
            zero_q1[quadrant->gt (ii)] = 0.;
            zeta0[quadrant->gt (ii)] = 1.0;
            zeta1[quadrant->gt (ii)] = 1.0;
            g1[quadrant->gt (ii)] = 0.;

            sold[ord[0](quadrant->gt (ii))] = 0.0;
            sold[ord[1](quadrant->gt (ii))] = 0.0;
            sold[ord[2](quadrant->gt (ii))] = 0.0;
            sold[ord[3](quadrant->gt (ii))] = 0.0;
            sold[ord[4](quadrant->gt (ii))] = 0.0;

            sol[ord[0](quadrant->gt (ii))] = 0.0;
            sol[ord[1](quadrant->gt (ii))] = 0.0;
            sol[ord[2](quadrant->gt (ii))] = 0.0;
            sol[ord[3](quadrant->gt (ii))] = 0.0;
            sol[ord[4](quadrant->gt (ii))] = 0.0;
          }
          else
            for (int jj = 0; jj < 2; ++jj)
              {
                zero_q1[quadrant->gparent (jj, ii)] += 0.;
                zeta0[quadrant->gparent (jj, ii)] += 0.;
                zeta1[quadrant->gparent (jj, ii)] += 0.;
                g1[quadrant->gparent (jj, ii)] += 0.;

                sold[ord[0](quadrant->gparent (jj, ii))] += 0.;
                sold[ord[1](quadrant->gparent (jj, ii))] += 0.;
                sold[ord[2](quadrant->gparent (jj, ii))] += 0.;
                sold[ord[3](quadrant->gparent (jj, ii))] += 0.;
                sold[ord[4](quadrant->gparent (jj, ii))] += 0.;
              }
        }
    }
//  tmsh.octbin_export_quadrant ("epsilon_file", epsilon);
//  tmsh.octbin_export_quadrant ("sigma_file", sigma);

  bim2a_solution_with_ghosts (tmsh, sold, replace_op, ord[0], false);
  bim2a_solution_with_ghosts (tmsh, sold, replace_op, ord[1], false);
  bim2a_solution_with_ghosts (tmsh, sold, replace_op, ord[2], false);
  bim2a_solution_with_ghosts (tmsh, sold, replace_op, ord[3], false);
  bim2a_solution_with_ghosts (tmsh, sold, replace_op, ord[4]);

  zero_q1.assemble (replace_op);
  zeta0.assemble (replace_op);
  zeta1.assemble (replace_op);
  g1.assemble (replace_op);                      

  // Save inital conditions
  sprintf(filename, "model_1_rho_0000");
  tmsh.octbin_export (filename, sold, ord[0]);
  sprintf(filename, "model_1_phi_0000");
  tmsh.octbin_export (filename, sold, ord[1]);

  sprintf(filename, "model_1_p1_0000");
  tmsh.octbin_export (filename, sold, ord[2]);
  sprintf(filename, "model_1_p2_0000");
  tmsh.octbin_export (filename, sold, ord[2]);
  sprintf(filename, "model_1_p3_0000");
  tmsh.octbin_export (filename, sold, ord[2]);

  int count = 0;

  //Choosing the indices for the nodes corresponding to the vales of rho of interest 
  /*
  if (rank == 0) {
    rho_out[0]=std::vector<double>(N_rhos+1,0.0);
    rho_idx = find_idx(tmsh,points,tols,N_rhos);
  }
  */

  // Time cycle
  double time = 0.0;
  double time_in_step = 0.0;
  double dt;
  double eps = 1.0e-10;
  double err_max;
  
  std::ofstream error_file;
  error_file.setf(std::ios::fixed);
  if (rank == 0 && save_error)
    error_file.open("error.txt");
  q1_vec sold1 = sold, sold2 = sold, sol1 = sol, sol2 = sol;
  while (time < T) {
    if (rank == 0)
      std::cout << "____________ COMPUTING FOR TIME = " << time + DT << " ____________" << std::endl;
    dt = DT;
    time_in_step = 0.0;
    while (true) {
      sold1 = sold; sold2 = sold; sol1 = sol; sol2 = sol;
      time_step(rank, time + time_in_step, dt,
                  ord, tmsh, lin_solver, A,
                  xa, ir, jc, epsilon, sigma,
                  zero_std_vect, zero_q1, delta1,delta0,
                  reaction_term_p1, reaction_term_p2, reaction_term_p3,
                  diffusion_term_p1, diffusion_term_p2, diffusion_term_p3,
                  zeta0, zeta1, f1, f0, g1, g0, gp1, gp2, gp3, sold1, sol1);
      time_step(rank, time + time_in_step, dt/2,
                  ord, tmsh, lin_solver, A,
                  xa, ir, jc, epsilon, sigma,
                  zero_std_vect, zero_q1, delta1,delta0,
                  reaction_term_p1, reaction_term_p2, reaction_term_p3,
                  diffusion_term_p1, diffusion_term_p2, diffusion_term_p3,
                  zeta0, zeta1, f1, f0, g1, g0, gp1, gp2, gp3, sold2, sol2);
      time_step(rank, time + time_in_step + dt/2, dt/2,
                  ord, tmsh, lin_solver, A,
                  xa, ir, jc, epsilon, sigma,
                  zero_std_vect, zero_q1, delta1,delta0,
                  reaction_term_p1, reaction_term_p2, reaction_term_p3,
                  diffusion_term_p1, diffusion_term_p2, diffusion_term_p3,
                  zeta0, zeta1, f1, f0, g1, g0, gp1, gp2, gp3, sold2, sol2);
      err_max = 0.;
      for (size_t i = 0; i < sold.local_size(); i++)
        err_max = std::max(err_max, std::fabs(sold1.get_owned_data()[i] - sold2.get_owned_data()[i]));
      
      MPI_Allreduce(MPI_IN_PLACE, &err_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

      if (rank == 0)
        std::cout << "error/toll = " << err_max/toll << std::endl;
      
      if (err_max < toll) {
        time_in_step += dt;
        sold = std::move(sold2);
        if (rank == 0 && save_error)
          error_file << "time = " << std::setw(12) << std::setprecision(5) << time + time_in_step << "    max err = " << std::setprecision(7) << err_max << std::endl;
        if (time_in_step > DT - eps) {
          time += DT;
          // Save solution
          if (save_sol == true) {
            ++count;
            sprintf(filename, "model_1_rho_%4.4d",count);
            tmsh.octbin_export (filename, sold, ord[0]);
            sprintf(filename, "model_1_phi_%4.4d",count);
            tmsh.octbin_export (filename, sold, ord[1]);
            sprintf(filename, "model_1_p1_%4.4d", count);
            tmsh.octbin_export (filename,sold, ord[2]);
            sprintf(filename, "model_1_p2_%4.4d", count);
            tmsh.octbin_export (filename,sold, ord[2]);
            sprintf(filename, "model_1_p3_%4.4d", count);
            tmsh.octbin_export (filename,sold, ord[2]);
          }
          break;
        }
        // update dt
        dt = std::min(DT-time_in_step, dt*std::pow(err_max/toll,-0.5)*0.9);
        
      }
      else
        dt/=2;
    }
  }
  error_file.close();
  /*
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
  */
  // Close MPI and print report
  MPI_Barrier (MPI_COMM_WORLD);

  // Clean linear solver
  lin_solver->cleanup ();
  
  MPI_Finalize (); 

  return 0;
}
