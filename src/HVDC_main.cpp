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
#include <vector>
#include <unordered_set>
#include <functional>

#include <bim_distributed_vector.h>
#include <bim_sparse_distributed.h>
#include <bim_timing.h>
#include <mumps_class.h>
#include <quad_operators_3d.h>

#include <nlohmann/json.hpp>


#include "plugins/test_factory.h"
#include "plugins/voltage_factory.h"
#include <dlfcn.h>
#include <filesystem>




extern double epsilon_0;




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

void json_export (std::ifstream &is, std::ofstream &os) {
  json J;
  std::vector<std::string> variable_names;
  size_t num_var = 0;
  std::string line;
  std::getline(is, line);
  std::stringstream header(line);
  while (!header.eof()){
    variable_names.push_back("");
    header >> variable_names[num_var];
    J[variable_names[num_var]] = std::vector<double>();
    num_var++;
  }

  double num;
  size_t count;
  while (!is.eof()){
    for (size_t i = 0; i < num_var; ++i){
      if (is >> num)
        J[variable_names[i]].push_back(num);
    }
  }
  os << std::setw(4) << J;
  return;
}

using q1_vec_ = q1_vec<distributed_vector>;
template <size_t N_eqs>
void time_step (const int rank, const double time, const double DELTAT,
                std::unique_ptr<tests::generic_test> const &test,
                std::unique_ptr<voltages::generic_voltage> const &voltage,
                const std::array<ordering,N_eqs> &ord,
                tmesh_3d &tmsh, mumps *lin_solver,
                distributed_sparse_matrix &A,
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
    dirichlet_bcs3 bcs0, bcs1;
    bcs0.push_back (std::make_tuple (0, 4, [&voltage, time, DELTAT](double x, double y, double z){return voltage->V_in_time(4, time+DELTAT, x, y, z);})); //bottom
    bcs0.push_back (std::make_tuple (0, 5, [&voltage, time, DELTAT](double x, double y, double z){return voltage->V_in_time(5, time+DELTAT, x, y, z);})); //top

    // Print curent time
    if (rank==0)
      std::cout << "TIME= " << time + DELTAT << std::endl;

    // Reset containers
    A.reset ();
    sol.get_owned_data().assign(sol.get_owned_data ().size (), 0.0);
    sol.assemble (replace_op);

    // Initialize non constant (in time) parameters
    for (auto quadrant = tmsh.begin_quadrant_sweep ();
    quadrant != tmsh.end_quadrant_sweep ();
    ++quadrant)
    {
      double xx{quadrant->centroid(0)}, yy{quadrant->centroid(1)}, zz{quadrant->centroid(2)};  

      reaction_term_p1[quadrant->get_forest_quad_idx ()] =
        1 + DELTAT / test->tau_p1_fun(xx,yy,zz);
      reaction_term_p2[quadrant->get_forest_quad_idx ()] =
        1 + DELTAT / test->tau_p2_fun(xx,yy,zz);
      reaction_term_p3[quadrant->get_forest_quad_idx ()] =
        1 + DELTAT / test->tau_p3_fun(xx,yy,zz);
      diffusion_term_p1[quadrant->get_forest_quad_idx ()] =
        - DELTAT / test->tau_p1_fun(xx,yy,zz) * epsilon_0 * test->csi_1_fun(xx,yy,zz);
      diffusion_term_p2[quadrant->get_forest_quad_idx ()] =
        - DELTAT / test->tau_p2_fun(xx,yy,zz) * epsilon_0 * test->csi_2_fun(xx,yy,zz);
      diffusion_term_p3[quadrant->get_forest_quad_idx ()] =
        - DELTAT / test->tau_p3_fun(xx,yy,zz) * epsilon_0 * test->csi_3_fun(xx,yy,zz);
      sigma[quadrant->get_forest_quad_idx ()] = test->sigma_fun(xx,yy,zz,DELTAT);
      for (int ii = 0; ii < 8; ++ii) {
        if (! quadrant->is_hanging (ii)){
          g0[quadrant->gt (ii)] = sold[ord[0](quadrant->gt (ii))];
          gp1[quadrant->gt (ii)] = sold[ord[2](quadrant->gt (ii))];
          gp2[quadrant->gt (ii)] = sold[ord[3](quadrant->gt (ii))];
          gp3[quadrant->gt (ii)] = sold[ord[4](quadrant->gt (ii))];
        }
        else
          for (int jj = 0; jj < quadrant->num_parents (ii); ++jj){
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
    
    // Advection_diffusion
    bim3a_advection_diffusion (tmsh, sigma, zero_q1, A, true, ord[0], ord[1]);
    bim3a_advection_diffusion (tmsh, epsilon, zero_q1, A, true, ord[1], ord[1]);
    bim3a_advection_diffusion (tmsh, diffusion_term_p1, zero_q1, A, true, ord[2], ord[1]);
    bim3a_advection_diffusion (tmsh, diffusion_term_p2, zero_q1, A, true, ord[3], ord[1]);
    bim3a_advection_diffusion (tmsh, diffusion_term_p3, zero_q1, A, true, ord[4], ord[1]);
    
    // Reaction
    bim3a_reaction (tmsh, delta0, zeta0, A, ord[0], ord[0]);
    bim3a_reaction (tmsh, delta1, zeta1, A, ord[1], ord[0]);
    bim3a_reaction (tmsh, delta0, zeta0, A, ord[1], ord[2]);
    bim3a_reaction (tmsh, delta0, zeta0, A, ord[1], ord[3]);
    bim3a_reaction (tmsh, delta0, zeta0, A, ord[1], ord[4]);
    bim3a_reaction (tmsh, reaction_term_p1, zeta1, A, ord[2], ord[2]);
    bim3a_reaction (tmsh, reaction_term_p2, zeta1, A, ord[3], ord[3]);
    bim3a_reaction (tmsh, reaction_term_p3, zeta1, A, ord[4], ord[4]);

    // Rhs
    bim3a_rhs (tmsh, f0, g0, sol, ord[0]);
    bim3a_rhs (tmsh, f1, g1, sol, ord[1]);
    bim3a_rhs (tmsh, f0, gp1, sol, ord[2]);
    bim3a_rhs (tmsh, f0, gp2, sol, ord[3]);
    bim3a_rhs (tmsh, f0, gp3, sol, ord[4]);

    // Boundary conditions
    bim3a_dirichlet_bc (tmsh, bcs0, A, sol, ord[0], ord[1], false);

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
  // alias definition
  using q1_vec = q1_vec<distributed_vector>;
  using json = nlohmann::json;

  // parsing data file
  std::ifstream data_file("data.json");
  json data = json::parse(data_file);
  data_file.close();

  // Initialize MPI
  MPI_Init (&argc, &argv);
  int rank, size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);
  
  std::string output_folder = data["output_location"];
  // select test to run and voltage function on contacts
  std::vector<std::string> test_name = data["test_to_run"];
  for (auto test_iter = test_name.cbegin(); test_iter != test_name.cend(); ++test_iter) {
  std::string vol_name = data[*test_iter]["algorithm"]["voltage_name"];
  
  // opening dynamic libraries
  std::string test_plugin = data[*test_iter]["physics_grid"]["physics_grid_plugin"];
  std::string voltage_plugin = data[*test_iter]["algorithm"]["voltage_plugin"];
  void *dl_test_p = dlopen(test_plugin.c_str(), RTLD_NOW);
  void *dl_voltage_p = dlopen(voltage_plugin.c_str(), RTLD_NOW);
  
  // importing some problem params
  double T = data[*test_iter]["algorithm"]["T"];
  epsilon_0 = data[*test_iter]["physics_grid"]["epsilon_0"];
  double Time = 0.; int count = 0;
  bool start_from_solution = data[*test_iter]["algorithm"]["start_from_solution"];
  int save_every_n_steps;
  std::string temp_solution_file_name;
  if (start_from_solution || data[*test_iter]["algorithm"]["save_temp_solution"])
    temp_solution_file_name = data[*test_iter]["algorithm"]["temp_sol"]["file"];
  bool save_temp_solution = data[*test_iter]["algorithm"]["save_temp_solution"];
  if (save_temp_solution) {
    save_every_n_steps = data[*test_iter]["algorithm"]["temp_sol"]["save_every_n_steps"];
  }
  if (start_from_solution) {
    Time = data[*test_iter]["algorithm"]["temp_sol"]["time"];
    count = data[*test_iter]["algorithm"]["temp_sol"]["count"];
  }
  
  double dt = data[*test_iter]["algorithm"]["initial_dt_for_adaptive_time_step"];
  double tol = data[*test_iter]["algorithm"]["tol_of_adaptive_time_step"];
  
  // set output preferences
  double DT = data[*test_iter]["options"]["print_solution_every_n_seconds"];
  bool save_sol = data[*test_iter]["options"]["save_sol"];
  bool save_error_and_comp_time = data[*test_iter]["options"]["save_error_and_comp_time"];
  bool save_currents = data[*test_iter]["options"]["save_currents"];
  bool save_displ_current = data[*test_iter]["options"]["save_displ_current"];
  bool compute_2_contacts = data[*test_iter]["options"]["compute_2_contacts"];

  // import from factories:
  // test
  auto which_test = tests::T_factory.find(*test_iter);
  std::unique_ptr<tests::generic_test> test = (which_test->second)();
  test->import_params(data);
  // voltage
  auto which_vol = voltages::V_factory.find(vol_name);
  std::unique_ptr<voltages::generic_voltage> voltage = (which_vol->second)();
  voltage->import_params(*test_iter, data);
  

  /*
  Manegement of solutions ordering:   Equation ordering:
  ord[0]-> rho                        ord[0]-> continuity equation
  ord[1]-> phi                        ord[1]-> diffusion-reaction equation
  ord[2]-> p1                         ord[2]-> p1 equation
  ord[3]-> p2                         ord[3]-> p2 equation
  ord[4]-> p3                         ord[4]-> p3 equation
  */
  
  constexpr size_t N_eqs= 5;
  constexpr size_t N_polcur = 3;
  const std::array<ordering,N_eqs> ord(makeorder<N_eqs>());
  const std::array<ordering,N_polcur+1> ord_c(makeorder<N_polcur+1>());
  const std::array<ordering,2> ord_displ_curr(makeorder<2>());


  // Generate the mesh in 3d
  tmesh_3d tmsh;
  tmsh.read_connectivity (simple_conn_p, simple_conn_num_vertices,
                          simple_conn_t, simple_conn_num_trees);

  // Uniform refinement
  int recursive = 1;
  tmsh.set_refine_marker ([&test](tmesh_3d::quadrant_iterator q) {return test->uniform_refinement(q);});
  tmsh.refine (recursive);

  // According to the test we may refine the grid or just leave it uniform
  if (test->extra_refinement) 
  {
    tmsh.set_refine_marker([&test](tmesh_3d::quadrant_iterator q) {return test->refinement(q);});
    tmsh.refine (recursive);

    tmsh.set_coarsen_marker([&test](tmesh_3d::quadrant_iterator q) {return test->coarsening(q);});
    tmsh.coarsen(recursive);
  }

  tmesh_3d::idx_t gn_nodes = tmsh.num_global_nodes ();
  tmesh_3d::idx_t ln_nodes = tmsh.num_owned_nodes ();
  tmesh_3d::idx_t ln_elements = tmsh.num_local_quadrants (); 

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

  // Declare matrix to compute currents
  distributed_sparse_matrix B;
  B.set_ranges (ln_nodes *2);

  // Buffer for export filename
  char filename[255]="";

  //Output rho vector
  //size_t N_timesteps = (size_t) (ceil(T/DELTAT)+1);
  //std::vector<std::vector<double>> rho_out(N_timesteps,std::vector<double>(N_rhos+1));

  // Compute coefficients

  // Diffusion
  std::vector<double> epsilon (ln_elements, 0.);
  std::vector<double> sigma (ln_elements, 0.);
  std::vector<double> zero_std_vect(ln_elements, 0.);
  std::vector<double> diffusion_term_p1 (ln_elements,0.);
  std::vector<double> diffusion_term_p2 (ln_elements,0.);
  std::vector<double> diffusion_term_p3 (ln_elements,0.);
  q1_vec zero_q1 (ln_nodes);

  // Reaction
  std::vector<double> delta1 (ln_elements, 0.);
  std::vector<double> delta0 (ln_elements, 0.);
  std::vector<double> reaction_term_p1 (ln_elements,0.);
  std::vector<double> reaction_term_p2 (ln_elements,0.);
  std::vector<double> reaction_term_p3 (ln_elements,0.);
  q1_vec zeta0 (ln_nodes);
  q1_vec zeta1 (ln_nodes);

  // Rhs
  std::vector<double> f1 (ln_elements, 0.);
  std::vector<double> f0 (ln_elements, 0.);
  std::vector<double> sigmaB (ln_elements, 0.);
  q1_vec g1 (ln_nodes);
  q1_vec g0 (ln_nodes);
  q1_vec gp1 (ln_nodes);
  q1_vec gp2 (ln_nodes);
  q1_vec gp3 (ln_nodes);

  q1_vec Bsol1 (ln_nodes * 2);
  q1_vec Bsol2 (ln_nodes * 2);
  q1_vec Ivec1 (ln_nodes * 2), Ivec2 (ln_nodes *2);
  std::unordered_set<size_t> Ivec_index1{};
  std::unordered_set<size_t> Ivec_index2{};
  
  // Setup streamings
  std::ofstream error_file, currents_file, I_displ_file;


  // Data to print
  std::array<double,4> charges;

  double I_c;
  double I_displ1, I_displ2, I_d1_c1, I_d2_c1, I_d1_c2, I_d2_c2;
  q1_vec Jx_vec(ln_nodes);
  q1_vec charges_vec(ln_nodes * (N_polcur+1));
  Jx_vec.get_owned_data().assign(Jx_vec.get_owned_data().size(),0.);
  charges_vec.get_owned_data().assign(charges_vec.get_owned_data().size(),0.);

  // Initialize constant (in time) parameters and initial data
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep ();
       ++quadrant) {
    double xx{quadrant->centroid(0)}, yy{quadrant->centroid(1)}, zz{quadrant->centroid(2)};  

    epsilon[quadrant->get_forest_quad_idx ()] = test->epsilon_fun(xx,yy,zz);
    delta1[quadrant->get_forest_quad_idx ()] = -1.0;
    delta0[quadrant->get_forest_quad_idx ()] = 1.0;
    f1[quadrant->get_forest_quad_idx ()] = 0.0;
    f0[quadrant->get_forest_quad_idx ()] = 1.0;
    sigmaB[quadrant->get_forest_quad_idx ()] = test->sigma_fun(xx,yy,zz,1.);

    for (int ii = 0; ii < 8; ++ii) {
      if (! quadrant->is_hanging (ii)) {
        zero_q1[quadrant->gt (ii)] = 0.;
        zeta0[quadrant->gt (ii)] = 1.0;
        zeta1[quadrant->gt (ii)] = 1.0;
        g1[quadrant->gt (ii)] = 0.;

        sol[ord[0](quadrant->gt (ii))] = 0.0;
        sol[ord[1](quadrant->gt (ii))] = 0.0;
        sol[ord[2](quadrant->gt (ii))] = 0.0;
        sol[ord[3](quadrant->gt (ii))] = 0.0;
        sol[ord[4](quadrant->gt (ii))] = 0.0;

        Jx_vec[quadrant->gt(ii)] = 0.;
        charges_vec[ord_c[0](quadrant->gt(ii))] = 0.;
        charges_vec[ord_c[1](quadrant->gt(ii))] = 0.;
        charges_vec[ord_c[2](quadrant->gt(ii))] = 0.;
        charges_vec[ord_c[3](quadrant->gt(ii))] = 0.;
      }
      else
        for (int jj = 0; jj < quadrant->num_parents (ii); ++jj) {
          zero_q1[quadrant->gparent (jj, ii)] += 0.;
          zeta0[quadrant->gparent (jj, ii)] += 0.;
          zeta1[quadrant->gparent (jj, ii)] += 0.;
          g1[quadrant->gparent (jj, ii)] += 0.;

          Jx_vec[quadrant->gparent (jj, ii)] += 0.;
          charges_vec[ord_c[0](quadrant->gparent (jj, ii))] += 0.;
          charges_vec[ord_c[1](quadrant->gparent (jj, ii))] += 0.;
          charges_vec[ord_c[2](quadrant->gparent (jj, ii))] += 0.;
          charges_vec[ord_c[3](quadrant->gparent (jj, ii))] += 0.;
        }
    }
  }
  sold.get_owned_data().assign(sold.local_size(),0.0);
  MPI_File temp_sol;
  if (start_from_solution) {
    MPI_File_open(MPI_COMM_WORLD, (temp_solution_file_name + "_" + std::to_string(count)).c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &temp_sol);
    MPI_File_seek(temp_sol, sold.get_range_start()*sizeof(double), MPI_SEEK_SET);
    MPI_File_read(temp_sol, sold.get_owned_data().data(), sold.local_size(), MPI_DOUBLE, MPI_STATUS_IGNORE);
    MPI_File_close(&temp_sol);
  }
  
  bim3a_solution_with_ghosts (tmsh, sold, replace_op, ord[0], false);
  bim3a_solution_with_ghosts (tmsh, sold, replace_op, ord[1], false);
  bim3a_solution_with_ghosts (tmsh, sold, replace_op, ord[2], false);
  bim3a_solution_with_ghosts (tmsh, sold, replace_op, ord[3], false);
  bim3a_solution_with_ghosts (tmsh, sold, replace_op, ord[4]);



  zero_q1.assemble (replace_op);
  zeta0.assemble (replace_op);
  zeta1.assemble (replace_op);
  g1.assemble (replace_op);

  // Save inital conditions
  if (!start_from_solution) {
    std::filesystem::create_directory(output_folder);
    std::filesystem::create_directory(output_folder + "/" + *test_iter);
    sprintf(filename, "%s/%s/model_1_rho_0000", output_folder.c_str(), test_iter->c_str());
    tmsh.octbin_export (filename, sold, ord[0]);
    sprintf(filename, "%s/%s/model_1_phi_0000", output_folder.c_str(), test_iter->c_str());
    tmsh.octbin_export (filename, sold, ord[1]);
    sprintf(filename, "%s/%s/model_1_p1_0000",  output_folder.c_str(), test_iter->c_str());
    tmsh.octbin_export (filename, sold, ord[2]);
    sprintf(filename, "%s/%s/model_1_p2_0000",  output_folder.c_str(), test_iter->c_str());
    tmsh.octbin_export (filename, sold, ord[3]);
    sprintf(filename, "%s/%s/model_1_p3_0000",  output_folder.c_str(), test_iter->c_str());
    tmsh.octbin_export (filename, sold, ord[4]);
  }


  // Lambda functions to use for currents computation (simple method)
  func3_quad Jx_mass = [&] (tmesh_3d::quadrant_iterator q, tmesh_3d::idx_t idx){
    return test->sigma_fun(q->centroid(0),q->centroid(1),q->centroid(2),1.)*
        (sold[ord[1](q->gt(4))] + sold[ord[1](q->gt(5))] + sold[ord[1](q->gt(6))] + sold[ord[1](q->gt(7))]
        -sold[ord[1](q->gt(0))] - sold[ord[1](q->gt(1))] - sold[ord[1](q->gt(2))] - sold[ord[1](q->gt(3))])/(q->p(2,4)-q->p(2,0))/4;
  };
  func3_quad free_charge_mass = [&] (tmesh_3d::quadrant_iterator q, tmesh_3d::idx_t idx){
    return (q->p(2,4)-q->p(2,0))*(sold[ord[0](q->gt(4))] + sold[ord[0](q->gt(5))] + sold[ord[0](q->gt(6))] + sold[ord[0](q->gt(7))])/8;
  };
  func3_quad p1_mass = [&] (tmesh_3d::quadrant_iterator q, tmesh_3d::idx_t idx){
    return (q->p(2,4)-q->p(2,0))*(sold[ord[2](q->gt(4))] + sold[ord[2](q->gt(5))] + sold[ord[2](q->gt(6))] + sold[ord[2](q->gt(7))])/8;
  };
  func3_quad p2_mass = [&] (tmesh_3d::quadrant_iterator q, tmesh_3d::idx_t idx){
    return (q->p(2,4)-q->p(2,0))*(sold[ord[3](q->gt(4))] + sold[ord[3](q->gt(5))] + sold[ord[3](q->gt(6))] + sold[ord[3](q->gt(7))])/8;
  };
  func3_quad p3_mass = [&] (tmesh_3d::quadrant_iterator q, tmesh_3d::idx_t idx){
    return (q->p(2,4)-q->p(2,0))*(sold[ord[4](q->gt(4))] + sold[ord[4](q->gt(5))] + sold[ord[4](q->gt(6))] + sold[ord[4](q->gt(7))])/8;
  };

  // Export test params
  if (rank == 0) {
    std::ofstream save_problem_data;
    save_problem_data.open(output_folder + "/" + *test_iter + "/" + "test.json");
    save_problem_data << std::setw(4) << data[*test_iter];
    save_problem_data.close();
  }
  // Print header of output files
  if (rank == 0 && save_error_and_comp_time) {
    if (!start_from_solution) {
      error_file.open(output_folder + "/" + *test_iter + "/" + "error_and_comp_time.txt");
      error_file << std::setw(20) << "time"
                 << std::setw(20) << "error/tol"
                 << std::setw(20) << "ts_comp_time"
                 << std::setw(20) << "total_time" << std::endl;
    }
    else
      error_file.open(output_folder + "/" + *test_iter + "/" + "error_and_comp_time.txt", std::fstream::app);
  
  }
  if (rank == 0 && save_currents) {
    if (!start_from_solution) {
      currents_file.open(output_folder + "/" + *test_iter + "/" + "currents_file.txt");
      if (! compute_2_contacts) {
        currents_file << std::setw(20) << "time"
                      << std::setw(20) << "I_c"
                      << std::setw(20) << "I_displ"
                      << std::setw(20) << "free_charge" 
                      << std::setw(20) << "P_inf_charge" 
                      << std::setw(20) << "P_1_charge" 
                      << std::setw(20) << "P_2_charge" 
                      << std::setw(20) << "P_3_charge" << std::endl;
      }
      else {
        currents_file << std::setw(20) << "time"
                      << std::setw(20) << "I_c"
                      << std::setw(20) << "I_displ1"
                      << std::setw(20) << "I_displ2"
                      << std::setw(20) << "free_charge" 
                      << std::setw(20) << "P_inf_charge" 
                      << std::setw(20) << "P_1_charge" 
                      << std::setw(20) << "P_2_charge" 
                      << std::setw(20) << "P_3_charge" << std::endl;
      }
    }
    else
      currents_file.open(output_folder + "/" + *test_iter + "/" + "currents_file.txt", std::fstream::app);
  }
  if (rank == 0 && save_displ_current) {
    if (!start_from_solution) {
      I_displ_file.open(output_folder + "/" + *test_iter + "/" + "I_displ_file.txt");
      if (!compute_2_contacts)
        I_displ_file  << std::setw(20) << "time"
                      << std::setw(20) << "I_displ" << std::endl;
      else
        I_displ_file  << std::setw(20) << "time"
                      << std::setw(20) << "I_displ1"
                      << std::setw(20) << "I_displ2" << std::endl;
    }
    else
      I_displ_file.open(output_folder + "/" + *test_iter + "/" + "I_displ_file.txt", std::fstream::app);

  }
  q1_vec sold1 = sold, sold2 = sold, sol1 = sol, sol2 = sol;

  // Store indexes of nodes on border where to estimate current
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
      quadrant != tmsh.end_quadrant_sweep (); ++quadrant)
    for (int ii = 0; ii < 8; ii++) {
      if (quadrant->e(ii) == 5)
        Ivec_index1.insert(ord_displ_curr[0](quadrant->gt (ii)));
      else if (compute_2_contacts || quadrant->e(ii) == 4)
        Ivec_index2.insert(ord_displ_curr[0](quadrant->gt (ii)));
    }
  
  // Time cycle
  double time_in_step = 0.0;
  double eps = 1.0e-10;
  double err_max;
  std::string last_saved_solution = "";
  double start_time, time0, time1;
  bool exit_loop;

  if (rank == 0)
    start_time = MPI_Wtime();
  while (Time < T - eps) {
    exit_loop = false;
    if (rank == 0)
      std::cout << "____________ COMPUTING FOR TIME = " << Time + DT << " ____________" << std::endl;
    time_in_step = 0.0;
    while (!exit_loop) {
      if(rank == 0) {
        std::cout << "dt = " << dt << std::endl;
        time0 = MPI_Wtime();
      }
      sold1 = sold; sold2 = sold; sol1 = sol; sol2 = sol;
      time_step<N_eqs>(rank, Time + time_in_step, dt, test, voltage,
                  ord, tmsh, lin_solver, A,
                  xa, ir, jc, epsilon, sigma,
                  zero_std_vect, zero_q1, delta1,delta0,
                  reaction_term_p1, reaction_term_p2, reaction_term_p3,
                  diffusion_term_p1, diffusion_term_p2, diffusion_term_p3,
                  zeta0, zeta1, f1, f0, g1, g0, gp1, gp2, gp3, sold1, sol1);
      time_step<N_eqs>(rank, Time + time_in_step, dt/2, test, voltage,
                  ord, tmsh, lin_solver, A,
                  xa, ir, jc, epsilon, sigma,
                  zero_std_vect, zero_q1, delta1,delta0,
                  reaction_term_p1, reaction_term_p2, reaction_term_p3,
                  diffusion_term_p1, diffusion_term_p2, diffusion_term_p3,
                  zeta0, zeta1, f1, f0, g1, g0, gp1, gp2, gp3, sold2, sol2);
      time_step<N_eqs>(rank, Time + time_in_step + dt/2, dt/2, test, voltage,
                  ord, tmsh, lin_solver, A,
                  xa, ir, jc, epsilon, sigma,
                  zero_std_vect, zero_q1, delta1,delta0,
                  reaction_term_p1, reaction_term_p2, reaction_term_p3,
                  diffusion_term_p1, diffusion_term_p2, diffusion_term_p3,
                  zeta0, zeta1, f1, f0, g1, g0, gp1, gp2, gp3, sold2, sol2);
      err_max = 0.;

      for (auto quadrant = tmsh.begin_quadrant_sweep ();
          quadrant != tmsh.end_quadrant_sweep (); ++quadrant)
        for (int ii = 0; ii < 8; ii++) {
            Bsol1[ord_displ_curr[0](quadrant->gt (ii))] =
              (sold1[ord[0](quadrant->gt (ii))] - sold[ord[0](quadrant->gt (ii))]) / dt;
            Bsol2[ord_displ_curr[0](quadrant->gt (ii))] =
              (sold2[ord[0](quadrant->gt (ii))] - sold[ord[0](quadrant->gt (ii))]) / dt;
            Bsol1[ord_displ_curr[1](quadrant->gt (ii))] = sold1[ord[1](quadrant->gt (ii))];
            Bsol2[ord_displ_curr[1](quadrant->gt (ii))] = sold2[ord[1](quadrant->gt (ii))];
        }

      Bsol1.assemble(replace_op);
      Bsol2.assemble(replace_op);
      
      B.reset();
      bim3a_reaction (tmsh, delta0, zeta0, B, ord_displ_curr[0], ord_displ_curr[0]);
      bim3a_advection_diffusion (tmsh, sigmaB, zero_q1, B, true, ord_displ_curr[0], ord_displ_curr[1]);
      Ivec1 = B*Bsol1;
      Ivec2 = B*Bsol2;

      I_d1_c1 = 0.; I_d2_c1 = 0.; I_d1_c2 = 0.; I_d2_c2 = 0.;
      for (auto it = Ivec_index1.cbegin(); it != Ivec_index1.cend(); it++) {
        I_d1_c1 += Ivec1[*it];
        I_d2_c1 += Ivec2[*it];
      }
      if (compute_2_contacts) {
        for (auto it = Ivec_index2.cbegin(); it != Ivec_index2.cend(); it++) {
          I_d1_c2 += Ivec1[*it];
          I_d2_c2 += Ivec2[*it];
        }
      }
      std::cout << "c1 rank " << rank << "    " << I_d1_c1 << "   " << I_d2_c1 << std::endl;
      std::cout << "c2 rank " << rank << "    " << I_d1_c2 << "   " << I_d2_c2 << std::endl;
      MPI_Allreduce(MPI_IN_PLACE, &I_d1_c1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &I_d2_c1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      if (compute_2_contacts) {
        MPI_Allreduce(MPI_IN_PLACE, &I_d1_c2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &I_d2_c2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      }

      I_displ1 = I_d2_c1;
      I_displ2 = I_d2_c2;
      err_max = std::fabs((I_d2_c1-I_d1_c1) / (std::fabs(I_d2_c1) + 1.0e-50));
      
      if (rank == 0)
        std::cout << "error/tol = " << err_max/tol << std::endl;
      
      if (err_max < tol) {
        time_in_step += dt;
        sold = std::move(sold2);
        if (rank == 0 && save_error_and_comp_time) {
          time1 = MPI_Wtime();
          error_file << std::setw(20) << std::setprecision(5) << Time + time_in_step
                     << std::setw(20) << std::setprecision(7) << err_max
                     << std::setw(20) << std::setprecision(7) << time1 - time0;
        }
        if (rank == 0 && save_displ_current) {
          if (!compute_2_contacts)
            I_displ_file  << std::setw(20) << Time + time_in_step
                          << std::setw(20) << I_displ1 << std::endl;
          else
            I_displ_file  << std::setw(20) << Time + time_in_step
                          << std::setw(20) << I_displ1
                          << std::setw(20) << I_displ2 << std::endl;
        }
        if (time_in_step > DT - eps) {
          Time += DT;
          ++count;
          // Save temp solution
          if (save_temp_solution && !(count % save_every_n_steps)) {
            MPI_File_open(MPI_COMM_WORLD, (temp_solution_file_name + "_" + std::to_string(count)).c_str(),
                          MPI_MODE_CREATE | MPI_MODE_WRONLY,
                          MPI_INFO_NULL, &temp_sol);
            MPI_File_seek(temp_sol, sold.get_range_start()*sizeof(double), MPI_SEEK_SET);
            MPI_File_write(temp_sol, sold.get_owned_data().data(),
                          sold.local_size(), MPI_DOUBLE, MPI_STATUS_IGNORE);
            MPI_File_close(&temp_sol);
            MPI_Barrier(MPI_COMM_WORLD);
            if (rank == 0)
              std::cout << "saved temp solution at time " + std::to_string(Time) 
                        << " and count " << count << std::endl;
            remove(last_saved_solution.c_str());
            last_saved_solution = temp_solution_file_name + "_" + std::to_string(count);
          }
          if (save_currents) {
            Jx_vec.get_owned_data().assign(Jx_vec.get_owned_data().size(), 0.);
            Jx_vec.assemble(replace_op);
            charges_vec.get_owned_data().assign(charges_vec.get_owned_data().size(), 0.);
            charges_vec.assemble(replace_op);

            bim3a_boundary_mass(tmsh, 0, 5, Jx_vec, Jx_mass);
            bim3a_boundary_mass(tmsh, 0, 5, charges_vec, free_charge_mass, ord_c[0]);
            bim3a_boundary_mass(tmsh, 0, 5, charges_vec, p1_mass, ord_c[1]);
            bim3a_boundary_mass(tmsh, 0, 5, charges_vec, p2_mass, ord_c[2]);
            bim3a_boundary_mass(tmsh, 0, 5, charges_vec, p3_mass, ord_c[3]);

            I_c = 0.;
            charges.fill(0.);

            for (size_t i = 0; i < Jx_vec.local_size(); i++)
              I_c += Jx_vec.get_owned_data()[i];
            for (size_t i = 0; i < charges_vec.local_size() / (N_polcur+1); i++) {
              charges[0] += charges_vec.get_owned_data()[i*(N_polcur+1)];
              charges[1] += charges_vec.get_owned_data()[i*(N_polcur+1)+1];
              charges[2] += charges_vec.get_owned_data()[i*(N_polcur+1)+2];
              charges[3] += charges_vec.get_owned_data()[i*(N_polcur+1)+3];
            }


            MPI_Allreduce(MPI_IN_PLACE, &I_c, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, charges.data(), 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            // I_p_inf = (E_flux - E_flux_old) / DT; E_flux_old = E_flux;
            
            if (!compute_2_contacts)
              currents_file << std::setw(20) << std::setprecision(5) << Time
                            << std::setw(20) << std::setprecision(5) << I_c
                            << std::setw(20) << std::setprecision(5) << I_displ1
                            << std::setw(20) << std::setprecision(5) << charges[0]
                            << std::setw(20) << std::setprecision(5) << - charges[0] - charges[1] - charges[2] - charges[3]
                            << std::setw(20) << std::setprecision(5) << charges[1]
                            << std::setw(20) << std::setprecision(5) << charges[2]
                            << std::setw(20) << std::setprecision(5) << charges[3] << std::endl;
            else
              currents_file << std::setw(20) << std::setprecision(5) << Time
                            << std::setw(20) << std::setprecision(5) << I_c
                            << std::setw(20) << std::setprecision(5) << I_displ1
                            << std::setw(20) << std::setprecision(5) << I_displ2
                            << std::setw(20) << std::setprecision(5) << charges[0]
                            << std::setw(20) << std::setprecision(5) << - charges[0] - charges[1] - charges[2] - charges[3]
                            << std::setw(20) << std::setprecision(5) << charges[1]
                            << std::setw(20) << std::setprecision(5) << charges[2]
                            << std::setw(20) << std::setprecision(5) << charges[3] << std::endl;
          }
          // Save solution
          if (save_sol == true) {
            sprintf(filename, "%s/%s/model_1_rho_%4.4d", output_folder.c_str(), test_iter->c_str(), count);
            tmsh.octbin_export (filename, sold, ord[0]);
            sprintf(filename, "%s/%s/model_1_phi_%4.4d", output_folder.c_str(), test_iter->c_str(), count);
            tmsh.octbin_export (filename, sold, ord[1]);
            sprintf(filename, "%s/%s/model_1_p1_%4.4d", output_folder.c_str(), test_iter->c_str(), count);
            tmsh.octbin_export (filename,sold, ord[2]);
            sprintf(filename, "%s/%s/model_1_p2_%4.4d", output_folder.c_str(), test_iter->c_str(), count);
            tmsh.octbin_export (filename,sold, ord[3]);
            sprintf(filename, "%s/%s/model_1_p3_%4.4d", output_folder.c_str(), test_iter->c_str(), count);
            tmsh.octbin_export (filename,sold, ord[4]);
          }
          dt = std::min(DT,dt*std::pow(err_max/tol,-0.5)*0.9);
          exit_loop = true;
        }
        // Update dt
        else {
          if (DT - time_in_step < dt*std::pow(err_max/tol,-0.5)*0.9)
            dt = DT-time_in_step;
          else
            dt = std::min((DT - time_in_step) / 2, dt*std::pow(err_max/tol,-0.5)*0.9);
        }
        if (rank == 0 && save_error_and_comp_time) {
            time1 = MPI_Wtime();
            error_file << std::setw(20) << std::setprecision(7) << time1 - start_time << std::endl;
          }
      }
      else
        dt /= 2;
    }
  }
  

  // Export in json format
  if (rank==0){
    error_file.close();
    currents_file.close();
    if(save_currents) {
      std::ifstream curr_file;
      std::ofstream curr_file_json;
      curr_file.open(output_folder + "/" + *test_iter + "/" + "currents_file.txt");
      curr_file_json.open(output_folder + "/" + *test_iter + "/" + "currents_file.json");
      json_export(curr_file, curr_file_json);
      curr_file.close();
      curr_file_json.close();
    }
    if (save_displ_current) {
      std::ifstream displ_curr_file;
      std::ofstream displ_curr_file_json;
      displ_curr_file.open(output_folder + "/" + *test_iter + "/" + "I_displ_file.txt");
      displ_curr_file_json.open(output_folder + "/" + *test_iter + "/" + "I_displ_file.json");
      json_export(displ_curr_file, displ_curr_file_json);
      displ_curr_file.close();
      displ_curr_file_json.close();
    }
    if (save_error_and_comp_time) {
      std::ifstream err_file;
      std::ofstream err_file_json;
      err_file.open(output_folder + "/" + *test_iter + "/" + "error_and_comp_time.txt");
      err_file_json.open(output_folder + "/" + *test_iter + "/" + "error_and_comp_time.json");
      json_export(err_file, err_file_json);
      err_file.close();
      err_file_json.close();
    }
  }
  dlclose(dl_test_p);
  dlclose(dl_voltage_p);
  // Close MPI and print report
  MPI_Barrier (MPI_COMM_WORLD);

  // Clean linear solver
  lin_solver->cleanup ();
  
  
  }
  MPI_Finalize (); 
  return 0;
}
