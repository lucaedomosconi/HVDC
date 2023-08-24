/*!
 * @file
 * @authors  Alessandro Lombardi, Luca Mosconi
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
#include <sstream>
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

//#include<test_2.h>
//#include<test_3.h>
//#include<test_4bis.h>
//#include<test_5.h>
#include "plugins/factory.h"
#include <dlfcn.h>
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

extern const int NUM_REFINEMENTS;
extern double T;
extern double tau;
extern double tau_p1, tau_p2, tau_p3;
extern bool save_sol;

		// Problem parameters
extern double epsilon_0;
constexpr bool fixed_sigma = false;
extern double tau_a;
extern double E_tr;

template<size_t N,class>
struct make_order_struct{};

template<size_t N, size_t... args>
struct make_order_struct<N, std::integer_sequence<size_t, args...>> {
  static auto fun() -> std::array<ordering, N> {
    return std::array<ordering, N> {dof_ordering<N, args>...};
  }
};

template<size_t N>
auto makeorder(){
  return make_order_struct<N,std::make_integer_sequence<size_t,N>>::fun();
}


void json_export(std::ifstream &is, std::ofstream &os) {
  json J;
  std::vector<std::string> variable_names;
  std::string prove;
  size_t num_var = 0;
  std::string line;
  std::getline(is, line);
  std::stringstream sstream(line);
  while (!sstream.eof()){
    variable_names.push_back("");
    sstream >> variable_names[num_var];
    J[variable_names[num_var]] = std::vector<double>();
    num_var++;
  }
  double num;
  size_t count;
  
  while(!is.eof()){
    for(size_t i = 0; i < num_var; ++i){
      if (is >> num)
        J[variable_names[i]].push_back(num);
    }
  }
  os << std::setw(4) << J;
  return;
}
using q1_vec_ = q1_vec<distributed_vector>;
template <size_t N_eqs>
void time_step(const int rank, const double time, const double DELTAT,
                std::unique_ptr<tests::generic_test> const &test,
                const std::array<ordering,N_eqs> &ord,
                tmesh &tmsh, mumps *lin_solver, distributed_sparse_matrix &A,
                std::vector<double> &xa, std::vector<int> &ir, std::vector<int> &jc,
                std::vector<double> &epsilon, std::vector<double> &sigma,
                std::vector<double> &zero_std_vect, q1_vec_ &zero_q1,
                std::vector<double> &delta1, std::vector<double> &delta0,
                std::vector<double> &reaction_term_p1,
                std::vector<double> &reaction_term_p2,
                std::vector<double> &reaction_term_p3,
                std::vector<double> &reaction_term_a,
                std::vector<double> &diffusion_term_p1,
                std::vector<double> &diffusion_term_p2,
                std::vector<double> &diffusion_term_p3,
                q1_vec_ &zeta0, q1_vec_ &zeta1,
                std::vector<double> &f1, std::vector<double> &f0,
                q1_vec_ &g1, q1_vec_ &g0,  q1_vec_ &gphi, q1_vec_ &gp1, q1_vec_ &gp2, q1_vec_ &gp3, q1_vec_ &ga1,
                gradient<q1_vec_> &grad_recovered, q1_vec_ &sold, q1_vec_ &sol) {

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
        - DELTAT / tau_p1 * epsilon_0 * test->csi_1_fun(xx,yy);
      diffusion_term_p2[quadrant->get_forest_quad_idx ()] =
        - DELTAT / tau_p2 * epsilon_0 * test->csi_2_fun(xx,yy);
      diffusion_term_p3[quadrant->get_forest_quad_idx ()] =
        - DELTAT / tau_p3 * epsilon_0 * test->csi_3_fun(xx,yy);
      if constexpr (fixed_sigma) 
        sigma[quadrant->get_forest_quad_idx ()] = test->sigma_fun(xx,yy,0.,0.,DELTAT);
      else {
        double grad_phi_module = std::pow(std::pow((sold[ord[1](quadrant->gt(1))]+sold[ord[1](quadrant->gt(3))]-sold[ord[1](quadrant->gt(0))]-sold[ord[1](quadrant->gt(2))])/(quadrant->p(0,1)-quadrant->p(0,0))/2,2)
                                        + std::pow((sold[ord[1](quadrant->gt(2))]+sold[ord[1](quadrant->gt(3))]-sold[ord[1](quadrant->gt(0))]-sold[ord[1](quadrant->gt(1))])/(quadrant->p(0,1)-quadrant->p(0,0))/2,2),0.5);
        double a_local = 0.25*(sold[ord[5](quadrant->gt(0))] + sold[ord[5](quadrant->gt(1))] + sold[ord[5](quadrant->gt(2))] + sold[ord[5](quadrant->gt(3))]);
        sigma[quadrant->get_forest_quad_idx ()] = test->sigma_fun(xx,yy,a_local,grad_phi_module,DELTAT);
        reaction_term_a[quadrant->get_forest_quad_idx ()] = 1 + DELTAT / tau_a;
      }
      if constexpr (!fixed_sigma)
        for (auto n_quadrant = quadrant->begin_neighbor_sweep();
            n_quadrant != quadrant->end_neighbor_sweep(); ++n_quadrant) {
          for (int ii = 0; ii < 4; ++ii) {
            if (!n_quadrant->is_hanging(ii))
              gphi[n_quadrant->gt (ii)] = sold[ord[1](n_quadrant->gt (ii))];
            else
              for (int jj = 0; jj < 2; ++jj)
                gphi[n_quadrant->gparent(jj, ii)] += 0;
          }
        }
      for (int ii = 0; ii < 4; ++ii) {
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
    }

    g0.assemble(replace_op);
    gp1.assemble(replace_op);
    gp2.assemble(replace_op);
    gp3.assemble(replace_op);
    
    if constexpr (!fixed_sigma) {
      gphi.assemble(replace_op);
      grad_recovered = bim2c_quadtree_pde_recovered_gradient(tmsh, gphi);
      for (auto quadrant = tmsh.begin_quadrant_sweep ();
            quadrant != tmsh.end_quadrant_sweep (); ++quadrant) {
        for (int ii = 0; ii < 4; ii++) {
          if (! quadrant-> is_hanging (ii))
            ga1[quadrant->gt (ii)] = sold[ord[5](quadrant->gt (ii))]
                  + DELTAT / tau_a * (std::pow(std::pow(std::get<0>(grad_recovered)[quadrant->gt (ii)],2)
                            + std::pow(std::get<1>(grad_recovered)[quadrant->gt (ii)],2),0.5) < E_tr);
          else
            for (int jj; jj < 2; ++jj)
              ga1[quadrant->gparent (jj, ii)] += 0.;
        }
      }
      ga1.assemble(replace_op);
    }

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
    if constexpr (!fixed_sigma) 
      bim2a_reaction (tmsh, reaction_term_a, zeta1, A, ord[5], ord[5]);
    //rhs
    bim2a_rhs (tmsh, f0, g0, sol, ord[0]);
    bim2a_rhs (tmsh, f1, g1, sol, ord[1]);
    bim2a_rhs (tmsh, f0, gp1, sol, ord[2]);
    bim2a_rhs (tmsh, f0, gp1, sol, ord[3]);
    bim2a_rhs (tmsh, f0, gp1, sol, ord[4]);
    if constexpr (!fixed_sigma) 
      bim2a_rhs (tmsh, f0, ga1, sol, ord[5]);
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

}





int
main (int argc, char **argv)
{

  using q1_vec = q1_vec<distributed_vector>;
  using json = nlohmann::json;

  std::ifstream data_file("data.json");
  json data = json::parse(data_file);
  void *dl_p = dlopen("libTests.so", RTLD_NOW);
  
  std::string test_to_run = data["test_to_run"];
  
  T = data[test_to_run]["T"];
  tau = data[test_to_run]["tau"];
  tau_p1 = data[test_to_run]["tau_p1"];
  tau_p2 = data[test_to_run]["tau_p2"];
  tau_p3 = data[test_to_run]["tau_p3"];
  save_sol = data[test_to_run]["save_sol"];
  
  epsilon_0 = data[test_to_run]["epsilon_0"];

  double DT = data[test_to_run]["DT"];
  double toll = data[test_to_run]["toll_of_adaptive_time_step"];
  bool save_error = data[test_to_run]["save_error"];
  bool save_currents = data[test_to_run]["save_currents"];

  
//  std::cout << dlerror() << std::endl; // uncomment this line only to debug!
  auto where = tests::factory.find(test_to_run);
  std::unique_ptr<tests::generic_test> const& test = (where->second)();
//  std::cout << test->works() << std::endl;
  test->import_params(data);
  data_file.close();
  /*
  Manegement of solutions ordering: ord[0]-> phi                 ord[1]->rho
  Equation ordering: ord[0]->diffusion-reaction equation         ord[1]->continuity equation 
  */
  
//  const std::array<ordering,N_eqs> ord{dof_ordering<N_eqs,0>,
//                                      dof_ordering<N_eqs,1>};
  
  constexpr size_t N_eqs = fixed_sigma ? 5 : 6;
  const std::array<ordering,N_eqs> ord(makeorder<N_eqs>());
  constexpr size_t N_polcur = 3;
  const std::array<ordering,3> ord_c(makeorder<3>());

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
  tmsh.set_refine_marker ([&test](tmesh::quadrant_iterator q) {return test->uniform_refinement(q);});
  tmsh.refine (recursive);

  //In test 1 we only have uniform refinement, in all other cases we perform additional refinement
  if (test->extra_refinement) 
  {
    tmsh.set_refine_marker([&test](tmesh::quadrant_iterator q) {return test->refinement(q);});
    tmsh.refine (recursive);

    tmsh.set_coarsen_marker([&test](tmesh::quadrant_iterator q) {return test->coarsening(q);});
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
  std::vector<double> reaction_term_a (ln_elements,0.);
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
  q1_vec ga1 (ln_nodes);
  q1_vec gphi (ln_nodes);
  gradient<q1_vec> grad_recovered = std::make_pair(q1_vec(ln_nodes), q1_vec(ln_nodes));

  // Initialize constant (in time) parameters and initial data
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep ();
       ++quadrant)
    {
      double xx{quadrant->centroid(0)}, yy{quadrant->centroid(1)};  

      epsilon[quadrant->get_forest_quad_idx ()] = test->epsilon_fun(xx,yy);
      
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
            if constexpr (!fixed_sigma) {
              sold[ord[5](quadrant->gt (ii))] = 0.0;
              std::get<0>(grad_recovered)[quadrant->gt (ii)] = 0.0;
              std::get<1>(grad_recovered)[quadrant->gt (ii)] = 0.0;
            }

            sol[ord[0](quadrant->gt (ii))] = 0.0;
            sol[ord[1](quadrant->gt (ii))] = 0.0;
            sol[ord[2](quadrant->gt (ii))] = 0.0;
            sol[ord[3](quadrant->gt (ii))] = 0.0;
            sol[ord[4](quadrant->gt (ii))] = 0.0;
            if constexpr (!fixed_sigma)
              sol[ord[5](quadrant->gt (ii))] = 0.0;
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
                if constexpr (!fixed_sigma) {
                  sold[ord[5](quadrant->gparent (jj, ii))] += 0.0;
                  std::get<0>(grad_recovered)[quadrant->gparent (jj, ii)] += 0.0;
                  std::get<1>(grad_recovered)[quadrant->gparent (jj, ii)] += 0.0;
                }
                sol[ord[0](quadrant->gparent (jj, ii))] += 0.0;
                sol[ord[1](quadrant->gparent (jj, ii))] += 0.0;
                sol[ord[2](quadrant->gparent (jj, ii))] += 0.0;
                sol[ord[3](quadrant->gparent (jj, ii))] += 0.0;
                sol[ord[4](quadrant->gparent (jj, ii))] += 0.0;
                if constexpr (!fixed_sigma) {
                  sol[ord[5](quadrant->gparent (jj, ii))] += 0.0;
                }
              }
        }
    }
//  tmsh.octbin_export_quadrant ("epsilon_file", epsilon);
//  tmsh.octbin_export_quadrant ("sigma_file", sigma);

  bim2a_solution_with_ghosts (tmsh, sold, replace_op, ord[0], false);
  bim2a_solution_with_ghosts (tmsh, sold, replace_op, ord[1], false);
  bim2a_solution_with_ghosts (tmsh, sold, replace_op, ord[2], false);
  bim2a_solution_with_ghosts (tmsh, sold, replace_op, ord[3], false);
  if constexpr (fixed_sigma)
    bim2a_solution_with_ghosts (tmsh, sold, replace_op, ord[4]);
  else {
    bim2a_solution_with_ghosts (tmsh, sold, replace_op, ord[4], false);
    bim2a_solution_with_ghosts (tmsh, sold, replace_op, ord[5]);
  }

  zero_q1.assemble (replace_op);
  zeta0.assemble (replace_op);
  zeta1.assemble (replace_op);
  g1.assemble (replace_op);                      

  // Save inital conditions
  if constexpr (fixed_sigma) {
    sprintf(filename, "output/model_1_rho_0000");
    tmsh.octbin_export (filename, sold, ord[0]);
    sprintf(filename, "output/model_1_phi_0000");
    tmsh.octbin_export (filename, sold, ord[1]);

    sprintf(filename, "output/model_1_p1_0000");
    tmsh.octbin_export (filename, sold, ord[2]);
    sprintf(filename, "output/model_1_p2_0000");
    tmsh.octbin_export (filename, sold, ord[3]);
    sprintf(filename, "output/model_1_p3_0000");
    tmsh.octbin_export (filename, sold, ord[4]);
  }
  else {
    sprintf(filename, "output/model_2_rho_0000");
    tmsh.octbin_export (filename, sold, ord[0]);
    sprintf(filename, "output/model_2_phi_0000");
    tmsh.octbin_export (filename, sold, ord[1]);

    sprintf(filename, "output/model_2_p1_0000");
    tmsh.octbin_export (filename, sold, ord[2]);
    sprintf(filename, "output/model_2_p2_0000");
    tmsh.octbin_export (filename, sold, ord[3]);
    sprintf(filename, "output/model_2_p3_0000");
    tmsh.octbin_export (filename, sold, ord[4]);
    sprintf(filename, "output/model_2_a_0000");
    tmsh.octbin_export (filename, sold, ord[5]);
  }

  int count = 0;

  // Time cycle
  double time = 0.0;
  double time_in_step = 0.0;
  double dt;
  double eps = 1.0e-10;
  double err_max;
  
  double pol_charge = 0., pol_charge_old = 0.;
  double I_c;
  double I_p_inf;
  double I_p_k;
  double E_flux, E_flux_old = 0.;
  q1_vec Jx_vec(ln_nodes);
  q1_vec Ex_vec(ln_nodes);
  q1_vec p_vec(ln_nodes * N_polcur);
  Jx_vec.get_owned_data().assign(Jx_vec.get_owned_data().size(),0.);
  Ex_vec.get_owned_data().assign(Ex_vec.get_owned_data().size(),0.);
  p_vec.get_owned_data().assign(p_vec.get_owned_data().size(),0.);
  
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
    quadrant != tmsh.end_quadrant_sweep ();
    ++quadrant) {
    for (int ii = 0; ii < 4; ++ii) {
      if (! quadrant->is_hanging (ii)) {
        Jx_vec[quadrant->gt(ii)] = 0.;
        Ex_vec[quadrant->gt(ii)] = 0.;
        p_vec[ord_c[0](quadrant->gt(ii))] = 0.;
        p_vec[ord_c[1](quadrant->gt(ii))] = 0.;
        p_vec[ord_c[2](quadrant->gt(ii))] = 0.;
      }
      else
        for (int jj = 0; jj < 2; ++jj) {
          Jx_vec[quadrant->gparent (jj, ii)] += 0.;
          Ex_vec[quadrant->gparent (jj, ii)] += 0.;
          p_vec[ord_c[0](quadrant->gparent (jj, ii))] += 0.;
          p_vec[ord_c[1](quadrant->gparent (jj, ii))] += 0.;
          p_vec[ord_c[2](quadrant->gparent (jj, ii))] += 0.;
        }
    }
  }

  func_quad Jx_mass = [&] (tmesh::quadrant_iterator q, tmesh::idx_t idx){
    return test->sigma_fun(q->centroid(0),q->centroid(1),0.,0.,1.)*(sold[ord[1](q->gt(1))]+sold[ord[1](q->gt(3))]-sold[ord[1](q->gt(0))]-sold[ord[1](q->gt(2))])/(q->p(0,1)-q->p(0,0))/2;
  };
  func_quad Ex_mass = [&] (tmesh::quadrant_iterator q, tmesh::idx_t idx){
    return test->epsilon_fun(q->centroid(0),q->centroid(1))*(sold[ord[1](q->gt(1))]+sold[ord[1](q->gt(3))]-sold[ord[1](q->gt(0))]-sold[ord[1](q->gt(2))])/(q->p(0,1)-q->p(0,0))/2;
  };
  func_quad p1_mass = [&] (tmesh::quadrant_iterator q, tmesh::idx_t idx){
    return (q->p(0,1)-q->p(0,0))*(sold[ord[2](q->gt(1))] + sold[ord[2](q->gt(3))])/4;
  };
  func_quad p2_mass = [&] (tmesh::quadrant_iterator q, tmesh::idx_t idx){
    return (q->p(0,1)-q->p(0,0))*(sold[ord[3](q->gt(1))] + sold[ord[3](q->gt(3))])/4;
  };
  func_quad p3_mass = [&] (tmesh::quadrant_iterator q, tmesh::idx_t idx){
    return (q->p(0,1)-q->p(0,0))*(sold[ord[4](q->gt(1))] + sold[ord[4](q->gt(3))])/4;
  };
  std::ofstream error_file, currents_file;
  
  


  if (rank == 0 && save_error) {
    error_file.open("error.txt");
    error_file << std::setw(20) << "time" << std::setw(20) << "max_error" << std::endl;
  }
  if (rank == 0 && save_currents) {
    currents_file.open("currents_file.txt");
    currents_file << std::setw(20) << "time" << std::setw(20) << "I_c" << std::setw(20) << "I_p_inf" << std::setw(20) << "I_p_k" << std::endl;
  }

  q1_vec sold1 = sold, sold2 = sold, sol1 = sol, sol2 = sol;

  dt = DT;

  while (time < T) {
    if (rank == 0)
      std::cout << "____________ COMPUTING FOR TIME = " << time + DT << " ____________" << std::endl;
    time_in_step = 0.0;
    while (true) {
      if(rank == 0)
        std::cout << "dt = " << dt << std::endl;
      sold1 = sold; sold2 = sold; sol1 = sol; sol2 = sol;
      time_step<N_eqs>(rank, time + time_in_step, dt, test,
                  ord, tmsh, lin_solver, A,
                  xa, ir, jc, epsilon, sigma,
                  zero_std_vect, zero_q1, delta1,delta0,
                  reaction_term_p1, reaction_term_p2, reaction_term_p3, reaction_term_a,
                  diffusion_term_p1, diffusion_term_p2, diffusion_term_p3,
                  zeta0, zeta1, f1, f0, g1, g0, gphi, gp1, gp2, gp3, ga1, grad_recovered, sold1, sol1);
      time_step<N_eqs>(rank, time + time_in_step, dt/2, test,
                  ord, tmsh, lin_solver, A,
                  xa, ir, jc, epsilon, sigma,
                  zero_std_vect, zero_q1, delta1,delta0,
                  reaction_term_p1, reaction_term_p2, reaction_term_p3, reaction_term_a,
                  diffusion_term_p1, diffusion_term_p2, diffusion_term_p3,
                  zeta0, zeta1, f1, f0, g1, g0, gphi, gp1, gp2, gp3, ga1, grad_recovered, sold2, sol2);
      time_step<N_eqs>(rank, time + time_in_step + dt/2, dt/2, test,
                  ord, tmsh, lin_solver, A,
                  xa, ir, jc, epsilon, sigma,
                  zero_std_vect, zero_q1, delta1,delta0,
                  reaction_term_p1, reaction_term_p2, reaction_term_p3, reaction_term_a,
                  diffusion_term_p1, diffusion_term_p2, diffusion_term_p3,
                  zeta0, zeta1, f1, f0, g1, g0, gphi, gp1, gp2, gp3, ga1, grad_recovered, sold2, sol2);
      err_max = 0.;
      for (size_t i = 0; i < sold.local_size() / N_eqs; i++)
        err_max = std::max(err_max, std::fabs(sold1.get_owned_data()[i*N_eqs+1] - sold2.get_owned_data()[i*N_eqs+1]));
      
      MPI_Allreduce(MPI_IN_PLACE, &err_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

      if (rank == 0)
        std::cout << "error/toll = " << err_max/toll << std::endl;
      
      if (err_max < toll) {
        time_in_step += dt;
        sold = std::move(sold2);
        if (rank == 0 && save_error)
          error_file << std::setw(20) << std::setprecision(5) << time + time_in_step << std::setw(20) << std::setprecision(7) << err_max << std::endl;
        if (time_in_step > DT - eps) {
          time += DT;
          if (save_currents) {
            
            Jx_vec.get_owned_data().assign(Jx_vec.get_owned_data().size(), 0.);
            Jx_vec.assemble(replace_op);
            Ex_vec.get_owned_data().assign(Ex_vec.get_owned_data().size(), 0.);
            Ex_vec.assemble(replace_op);
            p_vec.get_owned_data().assign(p_vec.get_owned_data().size(), 0.);
            p_vec.assemble(replace_op);

            bim2a_boundary_mass(tmsh, 0, 1, Jx_vec, Jx_mass);
            bim2a_boundary_mass(tmsh, 0, 1, Ex_vec, Ex_mass);
            bim2a_boundary_mass(tmsh, 0, 1, p_vec, p1_mass, ord_c[0]);
            bim2a_boundary_mass(tmsh, 0, 1, p_vec, p2_mass, ord_c[1]);
            bim2a_boundary_mass(tmsh, 0, 1, p_vec, p3_mass, ord_c[2]);

            I_c = 0.;
            E_flux = 0.;
            pol_charge = 0.;

            for (size_t i = 0; i < Jx_vec.local_size(); i++)
              I_c += Jx_vec.get_owned_data()[i];
            for (size_t i = 0; i < Ex_vec.local_size(); i++)
              E_flux += Ex_vec.get_owned_data()[i];
            for (size_t i = 0; i < p_vec.local_size(); i++)
              pol_charge += p_vec.get_owned_data()[i];


            MPI_Allreduce(MPI_IN_PLACE, &I_c, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &E_flux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &pol_charge, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            I_p_inf = (E_flux - E_flux_old) / DT;
            I_p_k = (pol_charge - pol_charge_old) / DT;
            pol_charge_old = pol_charge;
            E_flux_old = E_flux;
            currents_file << std::setw(20) << std::setprecision(5) << time << std::setw(20) << std::setprecision(5) << I_c << std::setw(20) << std::setprecision(5) << I_p_inf << std::setw(20) << std::setprecision(5) << I_p_k << std::endl;
          }
          // Save solution
          if (save_sol == true) {
            ++count;
            if constexpr (fixed_sigma) {
              sprintf(filename, "output/model_1_rho_%4.4d",count);
              tmsh.octbin_export (filename, sold, ord[0]);
              sprintf(filename, "output/model_1_phi_%4.4d",count);
              tmsh.octbin_export (filename, sold, ord[1]);
              sprintf(filename, "output/model_1_p1_%4.4d", count);
              tmsh.octbin_export (filename,sold, ord[2]);
              sprintf(filename, "output/model_1_p2_%4.4d", count);
              tmsh.octbin_export (filename,sold, ord[3]);
              sprintf(filename, "output/model_1_p3_%4.4d", count);
              tmsh.octbin_export (filename,sold, ord[4]);
            }
            else {
              sprintf(filename, "output/model_2_rho_%4.4d",count);
              tmsh.octbin_export (filename, sold, ord[0]);
              sprintf(filename, "output/model_2_phi_%4.4d",count);
              tmsh.octbin_export (filename, sold, ord[1]);
              sprintf(filename, "output/model_2_p1_%4.4d", count);
              tmsh.octbin_export (filename,sold, ord[2]);
              sprintf(filename, "output/model_2_p2_%4.4d", count);
              tmsh.octbin_export (filename,sold, ord[3]);
              sprintf(filename, "output/model_2_p3_%4.4d", count);
              tmsh.octbin_export (filename,sold, ord[4]);
              sprintf(filename, "output/model_2_a_%4.4d", count);
              tmsh.octbin_export (filename,sold, ord[5]);
            }
          }
          break;
        }
        // update dt
        if (DT - time_in_step < dt*std::pow(err_max/toll,-0.5)*0.9)
          dt = DT-time_in_step;
        else
          dt = std::min((DT - time_in_step) / 2, dt*std::pow(err_max/toll,-0.5)*0.9);
      }
      else
        dt /= 2;
    }
  }
  if (rank==0){
    error_file.close();
    currents_file.close();
    if(save_currents) {
      std::ifstream curr_file;
      std::ofstream curr_file_json;
      curr_file.open("currents_file.txt");
      curr_file_json.open("currents_file.json");
      json_export(curr_file, curr_file_json);
      curr_file.close();
      curr_file_json.close();
    }
    if(save_error) {
      std::ifstream err_file;
      std::ofstream err_file_json;
      err_file.open("error.txt");
      err_file_json.open("error.json");
      json_export(err_file, err_file_json);
      err_file.close();
      err_file_json.close();
    }
  }
  // Close MPI and print report
  MPI_Barrier (MPI_COMM_WORLD);

  // Clean linear solver
  lin_solver->cleanup ();
  
  MPI_Finalize (); 

  return 0;
}
