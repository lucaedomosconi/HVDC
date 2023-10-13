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
#include <dlfcn.h>
#include <filesystem>

#include <bim_distributed_vector.h>
#include <bim_sparse_distributed.h>
#include <bim_timing.h>
#include <mumps_class.h>
#include <quad_operators_3d.h>

#include <nlohmann/json.hpp>
#include "export_json.h"
#include "datagen.h"

#include "plugins/generic_factory.h"




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



using q1_vector = q1_vec<distributed_vector>;
template <size_t N_eqs>
void time_step (const int rank, const double time, const double DELTAT,
                std::unique_ptr<tests::generic_test> const &test,
                std::unique_ptr<voltages::generic_voltage> const &voltage,
                const std::array<ordering,N_eqs> &ord,
                tmesh_3d &tmsh, mumps *lin_solver,
                distributed_sparse_matrix &A,
                std::vector<double> &xa, std::vector<int> &ir, std::vector<int> &jc,
                std::vector<double> &epsilon, std::vector<double> &sigma,
                std::vector<double> &zero_std_vect, q1_vector &zero_q1,
                std::vector<double> &delta1, std::vector<double> &delta0,
                std::vector<double> &reaction_term_p1,
                std::vector<double> &reaction_term_p2,
                std::vector<double> &reaction_term_p3,
                std::vector<double> &diffusion_term_p1,
                std::vector<double> &diffusion_term_p2,
                std::vector<double> &diffusion_term_p3,
                q1_vector &zeta0, q1_vector &zeta1,
                std::vector<double> &f1, std::vector<double> &f0,
                q1_vector &g1, q1_vector &g0, q1_vector &gp1, q1_vector &gp2, q1_vector &gp3,
                q1_vector &sold, q1_vector &sol) {

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
    q1_vector result = lin_solver->get_distributed_solution ();
    for (int idx = sold.get_range_start (); idx < sold.get_range_end (); ++idx)
      sold (idx) = result (idx);
    sold.assemble (replace_op);

}





int
main (int argc, char **argv)
{
  // Alias definition
  using json = nlohmann::json;
  using testfactory = Factory<tests::generic_test, std::function<std::unique_ptr<tests::generic_test>()>>;
  using voltagefactory = Factory<voltages::generic_voltage, std::function<std::unique_ptr<voltages::generic_voltage>()>>;

  // Initialize MPI
  MPI_Init (&argc, &argv);
  int rank=0, size=1;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  // Options from command line
  bool warnings_on = true;
  std::string datafilename = "data.json";
  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i], "--no-warnings") == 0)
      warnings_on = false;
    else if (strcmp(argv[i], "--export-json") == 0) {
      if ((argc - i == 1 || argv[++i][0] == '-') && rank == 0)
          std::clog << "missing or invalid argument for --output-json" << std::endl;
      else if (rank == 0) {
        std::string filename(argv[i]);
        if (filename.substr(filename.length()-4, filename.length()) == ".txt")
          filename = filename.substr(0,filename.length()-4);
        txt2json(std::string(argv[i]), filename+".json");
      }
      MPI_Finalize();
      return 0;
    }
    else if (strcmp(argv[i], "--generate-params") == 0) {
      if ((argc-i!=1 && argv[++i][0] != '-') && rank == 0)
        print_data(argv[i]);
      else if (rank == 0)
        print_data();
      MPI_Finalize();
      return 0;
    }
    else if (strcmp(argv[i], "-f") == 0) {
      if ((argc - i == 1 || argv[++i][0] == '-') && rank == 0)
        std::clog << "missing or invalid argument after -f" << std::endl;
      else
        datafilename = argv[i];
    }
    else {
      if (rank == 0)
        std::clog << "Invalid argument" << std::endl;
      MPI_Finalize();
      return 0;
    }
  }

  // Get factories address
  testfactory & T_factory = testfactory::Instance();
  voltagefactory & V_factory = voltagefactory::Instance();

  // Parsing data file
  std::ifstream data_file(datafilename);
  json data = json::parse(data_file);
  data_file.close();

  // Set output folder
  std::string output_folder;
  try {output_folder = std::string(data["output_location"]);}
  catch (...) {std::cerr << "Error: Unable to read object [output_location]" << std::endl; throw;}
  if (!output_folder.empty())
    output_folder = output_folder + "/";

  // Select test to run and voltage function on contacts
  std::vector<std::string> test_name;

  try{test_name = data["test_to_run"];}
  catch (...) {std::cerr << "Error: Unable to read object [test_to_run]" << std::endl; throw;}

  // Check overwritings
  if (warnings_on) {
    bool  start_from_solution,
          save_sol,
          compute_charges_on_border,
          save_displ_current,
          save_error_and_comp_time;
    int give_warning;
    std::vector<std::string> files;
    if (rank == 0) {
      for (auto test_iter = test_name.cbegin(); test_iter != test_name.cend(); ++test_iter) {
        try {start_from_solution = data[*test_iter]["algorithm"]["start_from_solution"];}
        catch (...) {std::cerr << "Error: Unable to read object ["+*test_iter+"][algorithm][start_from_solution]" << std::endl; throw;}
        try {save_sol = data[*test_iter]["options"]["save_sol"];}
        catch (...) {std::cerr << "Error: Unable to read object ["+*test_iter+"][options][save_sol]" << std::endl; throw;}
        try {compute_charges_on_border = data[*test_iter]["options"]["compute_charges_on_border"];}
        catch (...) {std::cerr << "Error: Unable to read object ["+*test_iter+"][options][compute_charges_on_border]" << std::endl; throw;}
        try {save_displ_current = data[*test_iter]["options"]["save_displ_current"];}
        catch (...) {std::cerr << "Error: Unable to read object ["+*test_iter+"][options][save_displ_current]" << std::endl; throw;}
        try {save_error_and_comp_time = data[*test_iter]["options"]["save_error_and_comp_time"];}
        catch (...) {std::cerr << "Error: Unable to read object ["+*test_iter+"][options][save_error_and_comp_time]" << std::endl; throw;}
        if (!start_from_solution) {
          if (save_error_and_comp_time  && std::filesystem::exists(output_folder + *test_iter + "/error_and_comp_time.txt")) {files.push_back(output_folder + *test_iter + "/error_and_comp_time.txt");}
          if (compute_charges_on_border && std::filesystem::exists(output_folder + *test_iter + "/charges_file.txt")) {files.push_back(output_folder + *test_iter + "/charges_file.txt");}
          if (save_displ_current        && std::filesystem::exists(output_folder + *test_iter + "/I_displ_file.txt")) {files.push_back(output_folder + *test_iter + "/I_displ_file.txt");}
          if (save_sol                  && std::filesystem::exists(output_folder + *test_iter + "/sol")) {files.push_back(output_folder + *test_iter + "/sol");}
        }
      }
      give_warning = files.empty() ? 0 : 1;
    }
    MPI_Bcast(&give_warning, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (give_warning) {
      char proceed = ' ';
      if (rank == 0) {
        std::clog << "The following files will be overwritten:" << std::endl;
        for (auto it = files.cbegin(); it != files.cend(); ++it)
          std::clog << *it << std::endl;
        std::clog << "Proceed anyway? ('y' or 'n')" << std::endl;
        while (proceed != 'y' && proceed != 'n') 
          std::cin >> proceed;
      }
      MPI_Bcast(&proceed, 1, MPI_CHAR, 0, MPI_COMM_WORLD);
      if (proceed == 'n') {
        MPI_Finalize();
        return 0;
      }
    }
  }
  // Iterate over the tests to run
  for (auto test_iter = test_name.cbegin(); test_iter != test_name.cend(); ++test_iter) {
  
  std::string vol_name;
  try{vol_name = data[*test_iter]["algorithm"]["voltage_name"];}
  catch(...) {std::cerr << "Error: Unable to read object: ["
                        << *test_iter
                        << "][algorithm][voltage_name]"
                        << std::endl; throw;}
  
  // Opening dynamic libraries
  std::string test_plugin;
  std::string voltage_plugin;
  try{test_plugin = data[*test_iter]["physics"]["physics_plugin"];}
  catch(...) {std::cerr << "Error: Unable to read object: ["
                        << *test_iter
                        << "][physics][physics_plugin]"
                        << std::endl; throw;}
  try{voltage_plugin = data[*test_iter]["algorithm"]["voltage_plugin"];}
  catch(...) {std::cerr << "Error: Unable to read object: ["
                        << *test_iter
                        << "][physics][voltage_plugin]"
                        << std::endl; throw;}
  void *dl_test_p = dlopen(test_plugin.c_str(), RTLD_NOW);
  void *dl_voltage_p = dlopen(voltage_plugin.c_str(), RTLD_NOW);
  
  // Importing some problem params
  // full time of the problem
  double T;
  // The physical constant
  epsilon_0;
  // Time of current step
  double Time;
  // max DT for adaptive time step, every DT the solution is printed
  double DT;
  // starting dt of adaptive time step
  double dt;
  // used to compute total computational time of a interrupted and resumed simulation
  double comp_time_of_previous_simuls;
  // tolerance of adaptiv time step
  double tol;
  // time step count
  int count;
  // how frequently save a temporary solution
  int save_every_n_steps;
  // list of options, reference to README.md
  bool start_from_solution;
  bool save_temp_solution;
  bool save_sol;
  bool save_error_and_comp_time;
  bool save_charges;
  bool save_displ_current;
  bool compute_2_contacts;
  // name of temporary solution file to resume the simulation with
  std::string temp_solution_file_name;

  Time = 0.; count = 0;
  try {T = data[*test_iter]["algorithm"]["T"];}
  catch (...) {std::cerr << "Error: Unable to read object ["+*test_iter+"][algorithm][T]" << std::endl; throw;}
  try {start_from_solution = data[*test_iter]["algorithm"]["start_from_solution"];}
  catch (...) {std::cerr << "Error: Unable to read object ["+*test_iter+"][algorithm][start_from_solution]" << std::endl; throw;}
  try {epsilon_0 = data[*test_iter]["physics"]["epsilon_0"];}
  catch (...) {std::cerr << "Error: Unable to read object ["+*test_iter+"][algorithm][epsilon_0]" << std::endl; throw;}
  try {save_temp_solution = data[*test_iter]["algorithm"]["save_temp_solution"];}
  catch (...) {std::cerr << "Error: Unable to read object ["+*test_iter+"][algorithm][save_temp_solution]" << std::endl; throw;}
  try {
    if (start_from_solution) {
      temp_solution_file_name = data[*test_iter]["algorithm"]["temp_sol"]["file_of_starting_sol"];
      if (!std::filesystem::exists(temp_solution_file_name)) {
        std::cerr << "Error: Temporary solution file \"" << temp_solution_file_name << " does not exist!" << std::endl;
        std::exit(1);
      }
    }
  }
  catch (...) {std::cerr << "Error: Unable to read object ["+*test_iter+"][algorithm][temp_sol][file_of_starting_sol]" << std::endl; throw;}
  
  try {
    if (save_temp_solution) {
      save_every_n_steps = data[*test_iter]["algorithm"]["temp_sol"]["save_every_n_steps"];
    }
  }
  catch(...) {std::cerr << "Error: Unable to read object ["+*test_iter+"][algorithm][temp_sol][sae_every_n_steps]" << std::endl; throw;}
  comp_time_of_previous_simuls = 0.0;

  try {dt = data[*test_iter]["algorithm"]["initial_dt_for_adaptive_time_step"];}
  catch (...) {std::cerr << "Error: Unable to read object ["+*test_iter+"][algorithm][initial_dt_for_adaptive_time_step]" << std::endl; throw;}
  try {tol = data[*test_iter]["algorithm"]["tol_of_adaptive_time_step"];}
  catch (...) {std::cerr << "Error: Unable to read object ["+*test_iter+"][algorithm][tol_of_adaptive_time_step]" << std::endl; throw;}
  // Set output preferences
  try {DT = data[*test_iter]["options"]["print_solution_every_n_seconds"];}
  catch (...) {std::cerr << "Error: Unable to read object ["+*test_iter+"][options][print_solution_every_n_seconds]" << std::endl; throw;}
  try {save_sol = data[*test_iter]["options"]["save_sol"];}
  catch (...) {std::cerr << "Error: Unable to read object ["+*test_iter+"][options][save_sol]" << std::endl; throw;}
  try {save_error_and_comp_time = data[*test_iter]["options"]["save_error_and_comp_time"];}
  catch (...) {std::cerr << "Error: Unable to read object ["+*test_iter+"][options][save_error_and_comp_time]" << std::endl; throw;}
  try {save_charges = data[*test_iter]["options"]["compute_charges_on_border"];}
  catch (...) {std::cerr << "Error: Unable to read object ["+*test_iter+"][options][compute_charges_on_border]" << std::endl; throw;}
  try {save_displ_current = data[*test_iter]["options"]["save_displ_current"];}
  catch (...) {std::cerr << "Error: Unable to read object ["+*test_iter+"][options][save_displ_current]" << std::endl; throw;}
  try {compute_2_contacts = data[*test_iter]["options"]["compute_2_contacts"];}
  catch (...) {std::cerr << "Error: Unable to read object ["+*test_iter+"][options][compute_2_contacts]" << std::endl; throw;}

  // Test
  std::unique_ptr<tests::generic_test> test = nullptr;
  try {test = T_factory.create(data[*test_iter]["physics"]["plugin_test_index"]);}
  catch (std::runtime_error const & e) {std::cerr << e.what() << std::endl; continue;}
  try {test->import_params(data[*test_iter]);}
  catch (std::runtime_error const & e) {std::cerr << "Error: Unable to read [" << *test_iter << "]" << e.what() <<std::endl;}
  catch (...) {std::cerr << "Check typos or missing elements among phisics plugin parameters" << *test_iter << std::endl; throw;}
  // Voltage
  std::unique_ptr<voltages::generic_voltage> voltage;
  try{voltage = V_factory.create(vol_name);}
  catch(std::runtime_error const & e) {std::cerr << e.what() << std::endl; continue;}
  try{voltage->import_params(data[*test_iter]["algorithm"]["voltage_plugin_params"]);}
  catch(std::runtime_error const & e) {std::cerr << "Error: Unable to read [" + *test_iter + "][algorithm][voltage_plugin_params]" << e.what() <<std::endl;}
  catch(...) {std::cerr << "Check typos or missing elements among voltage plugin parameters in " << *test_iter << std::endl; throw;}
  

  /*
  Manegement of solutions ordering:   Equation ordering:
  ord[0]-> rho                        ord[0]-> continuity equation
  ord[1]-> phi                        ord[1]-> diffusion-reaction equation
  ord[2]-> p1                         ord[2]-> p1 equation
  ord[3]-> p2                         ord[3]-> p2 equation
  ord[4]-> p3                         ord[4]-> p3 equation
  */
  
  // Setup number of equations, number of polarization currents and ordering arrays
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
  q1_vector sold (ln_nodes * N_eqs);
  sold.get_owned_data ().assign (sold.get_owned_data ().size (), 0.0);

  q1_vector sol (ln_nodes * N_eqs);
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

  // Compute coefficients

  // Diffusion
  std::vector<double> epsilon (ln_elements, 0.);
  std::vector<double> sigma (ln_elements, 0.);
  std::vector<double> zero_std_vect(ln_elements, 0.);
  std::vector<double> diffusion_term_p1 (ln_elements,0.);
  std::vector<double> diffusion_term_p2 (ln_elements,0.);
  std::vector<double> diffusion_term_p3 (ln_elements,0.);
  q1_vector zero_q1 (ln_nodes);

  // Reaction
  std::vector<double> delta1 (ln_elements, 0.);
  std::vector<double> delta0 (ln_elements, 0.);
  std::vector<double> reaction_term_p1 (ln_elements,0.);
  std::vector<double> reaction_term_p2 (ln_elements,0.);
  std::vector<double> reaction_term_p3 (ln_elements,0.);
  q1_vector zeta0 (ln_nodes);
  q1_vector zeta1 (ln_nodes);

  // Rhs
  std::vector<double> f1 (ln_elements, 0.);
  std::vector<double> f0 (ln_elements, 0.);
  std::vector<double> sigmaB (ln_elements, 0.);
  q1_vector g1 (ln_nodes);
  q1_vector g0 (ln_nodes);
  q1_vector gp1 (ln_nodes);
  q1_vector gp2 (ln_nodes);
  q1_vector gp3 (ln_nodes);

  // Variables for I_displ computation
  q1_vector Bsol1 (ln_nodes * 2);
  q1_vector Bsol2 (ln_nodes * 2);
  q1_vector Ivec1 (ln_nodes * 2), Ivec2 (ln_nodes *2);
  std::unordered_set<size_t> Ivec_index1{};
  std::unordered_set<size_t> Ivec_index2{};
  
  // Setup streamings
  std::ofstream error_file, charges_file, I_displ_file;
  const std::string error_file_name = output_folder + *test_iter + "/" + "error_and_comp_time.txt";
  const std::string charges_file_name = output_folder + *test_iter + "/" + "charges_file.txt";
  const std::string I_displ_file_name = output_folder + *test_iter + "/" + "I_displ_file.txt";

  // Data to print
  std::array<double,4> rho_pi_k;

  double I_c;
  double Ez_eps0;
  double I_displ1, I_displ2, I_d1_c1, I_d2_c1, I_d1_c2, I_d2_c2;
  q1_vector Jz_vec(ln_nodes), Ez_eps0_vec(ln_nodes);
  q1_vector rho_pi_k_vec(ln_nodes * (N_polcur+1));
  Jz_vec.get_owned_data().assign(Jz_vec.get_owned_data().size(),0.);
  Ez_eps0_vec.get_owned_data().assign(Ez_eps0_vec.get_owned_data().size(),0.);
  rho_pi_k_vec.get_owned_data().assign(rho_pi_k_vec.get_owned_data().size(),0.);

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

        Jz_vec[quadrant->gt(ii)] = 0.;
        Ez_eps0_vec[quadrant->gt(ii)] = 0.;
        rho_pi_k_vec[ord_c[0](quadrant->gt(ii))] = 0.;
        rho_pi_k_vec[ord_c[1](quadrant->gt(ii))] = 0.;
        rho_pi_k_vec[ord_c[2](quadrant->gt(ii))] = 0.;
        rho_pi_k_vec[ord_c[3](quadrant->gt(ii))] = 0.;
      }
      else
        for (int jj = 0; jj < quadrant->num_parents (ii); ++jj) {
          zero_q1[quadrant->gparent (jj, ii)] += 0.;
          zeta0[quadrant->gparent (jj, ii)] += 0.;
          zeta1[quadrant->gparent (jj, ii)] += 0.;
          g1[quadrant->gparent (jj, ii)] += 0.;

          Jz_vec[quadrant->gparent (jj, ii)] += 0.;
          Ez_eps0_vec[quadrant->gparent (jj, ii)] += 0.;
          rho_pi_k_vec[ord_c[0](quadrant->gparent (jj, ii))] += 0.;
          rho_pi_k_vec[ord_c[1](quadrant->gparent (jj, ii))] += 0.;
          rho_pi_k_vec[ord_c[2](quadrant->gparent (jj, ii))] += 0.;
          rho_pi_k_vec[ord_c[3](quadrant->gparent (jj, ii))] += 0.;
        }
    }
  }
  sold.get_owned_data().assign(sold.local_size(),0.0);
  
  // Read temp. solution
  MPI_File temp_sol;
  if (start_from_solution) {
    MPI_File_open(MPI_COMM_WORLD, temp_solution_file_name.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &temp_sol);
    if (rank == 0) {
      MPI_File_read_at(temp_sol, 0, &count, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
      MPI_File_read_at(temp_sol, sizeof(int), &Time, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
      MPI_File_read_at(temp_sol, sizeof(int)+sizeof(double), &comp_time_of_previous_simuls, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
      std::clog << "starting from solution_file \"" << temp_solution_file_name
                << "\"\nAt time = " << Time << ";"
                << "\ncount = " << count << std::endl;
    }
    MPI_File_read_at(temp_sol, sizeof(int)+(sold.get_range_start()+2)*sizeof(double), sold.get_owned_data().data(), sold.local_size(), MPI_DOUBLE, MPI_STATUS_IGNORE);
    MPI_File_close(&temp_sol);
    MPI_Bcast(&Time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&count, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }
  temp_solution_file_name = *test_iter + "_temp_sol";
  
  bim3a_solution_with_ghosts (tmsh, sold, replace_op, ord[0], false);
  bim3a_solution_with_ghosts (tmsh, sold, replace_op, ord[1], false);
  bim3a_solution_with_ghosts (tmsh, sold, replace_op, ord[2], false);
  bim3a_solution_with_ghosts (tmsh, sold, replace_op, ord[3], false);
  bim3a_solution_with_ghosts (tmsh, sold, replace_op, ord[4]);



  zero_q1.assemble (replace_op);
  zeta0.assemble (replace_op);
  zeta1.assemble (replace_op);
  g1.assemble (replace_op);

  // Create directories for output
  if (!output_folder.empty())
    std::filesystem::create_directory(output_folder);
  std::filesystem::create_directory(output_folder + *test_iter);
  std::filesystem::create_directory(output_folder + *test_iter + "/" + "sol");

  MPI_Barrier(MPI_COMM_WORLD);
  // Save inital conditions
  if (save_sol && !start_from_solution) {
    sprintf(filename, "%s%s/%s/model_1_rho_0000", output_folder.c_str(), test_iter->c_str(), "sol");
    tmsh.octbin_export (filename, sold, ord[0]);
    sprintf(filename, "%s%s/%s/model_1_phi_0000", output_folder.c_str(), test_iter->c_str(), "sol");
    tmsh.octbin_export (filename, sold, ord[1]);
    sprintf(filename, "%s%s/%s/model_1_p1_0000",  output_folder.c_str(), test_iter->c_str(), "sol");
    tmsh.octbin_export (filename, sold, ord[2]);
    sprintf(filename, "%s%s/%s/model_1_p2_0000",  output_folder.c_str(), test_iter->c_str(), "sol");
    tmsh.octbin_export (filename, sold, ord[3]);
    sprintf(filename, "%s%s/%s/model_1_p3_0000",  output_folder.c_str(), test_iter->c_str(), "sol");
    tmsh.octbin_export (filename, sold, ord[4]);
  }



  // Functions to use for currents computation (simple method)
  func3_quad Jz_mass = [&] (tmesh_3d::quadrant_iterator q, tmesh_3d::idx_t idx){
    return test->sigma_fun(q->centroid(0),q->centroid(1),q->centroid(2),1.)*
        (sold[ord[1](q->gt(4))] + sold[ord[1](q->gt(5))] + sold[ord[1](q->gt(6))] + sold[ord[1](q->gt(7))]
        -sold[ord[1](q->gt(0))] - sold[ord[1](q->gt(1))] - sold[ord[1](q->gt(2))] - sold[ord[1](q->gt(3))])/(q->p(2,4)-q->p(2,0))/4;
  };
  
  func3_quad Ez_eps0_mass = [&] (tmesh_3d::quadrant_iterator q, tmesh_3d::idx_t idx){
    return epsilon_0*
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
    save_problem_data.open(output_folder + *test_iter + "/" + "test.json");
    save_problem_data << std::setw(4) << data[*test_iter];
    save_problem_data.close();
  }

  // Print header of output files
  // ... for error_and_computation_file
  if (rank == 0 && save_error_and_comp_time) {
    if (!start_from_solution) {
      error_file.open(error_file_name);
      error_file << std::setw(20) << "time"
                 << std::setw(20) << "error"
                 << std::setw(20) << "ts_comp_time"
                 << std::setw(20) << "total_time" << std::endl;
    }
    else
      error_file.open(error_file_name, std::fstream::app);
  
  }

  // ... for charges_file
  if (rank == 0 && save_charges) {
    if (!start_from_solution) {
      charges_file.open(charges_file_name);
      if (! compute_2_contacts) {
        charges_file << std::setw(20) << "time"
                     << std::setw(20) << "I_c"
                     << std::setw(20) << "I_displ"
                     << std::setw(20) << "free_charge" 
                     << std::setw(20) << "P_inf_charge" 
                     << std::setw(20) << "P_1_charge" 
                     << std::setw(20) << "P_2_charge" 
                     << std::setw(20) << "P_3_charge" << std::endl;
      }
      else {
        charges_file << std::setw(20) << "time"
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
      charges_file.open(charges_file_name, std::fstream::app);
  }

  // For I_displ_file
  if (rank == 0 && save_displ_current) {
    if (!start_from_solution) {
      I_displ_file.open(I_displ_file_name);
      if (!compute_2_contacts)
        I_displ_file  << std::setw(20) << "time"
                      << std::setw(20) << "I_displ" << std::endl;
      else
        I_displ_file  << std::setw(20) << "time"
                      << std::setw(20) << "I_displ1"
                      << std::setw(20) << "I_displ2" << std::endl;
    }
    else
      I_displ_file.open(I_displ_file_name, std::fstream::app);

  }
  q1_vector sold1 = sold, sold2 = sold, sol1 = sol, sol2 = sol;

  // Store indexes of nodes on border where to estimate current
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
      quadrant != tmsh.end_quadrant_sweep (); ++quadrant)
    for (int ii = 0; ii < 8; ii++) {
      if (quadrant->e(ii) == 5) {
        Ivec_index1.insert(ord_displ_curr[0](quadrant->gt (ii)));
      }
      else if (compute_2_contacts && quadrant->e(ii) == 4) {
        Ivec_index2.insert(ord_displ_curr[0](quadrant->gt (ii)));
      }
    }
  
  // Time cycle
  double time_in_step = 0;
  double dt_start_big_step = dt;
  bool truncated_dt;
  double eps = 1.0e-8;
  double err_max;
  std::string last_saved_solution = "";
  double start_time, time0, time1;
  bool exit_loop;

  if (rank == 0)
    start_time = MPI_Wtime();
  while (Time < T - eps) {
    exit_loop = false;
    truncated_dt = 0;
    dt = dt_start_big_step;
    dt_start_big_step = 0.;
    if (rank == 0)
      std::cout << "____________ COMPUTING FOR TIME = " << Time + DT << " ____________" << std::endl;
    time_in_step = 0.0;
    while (!exit_loop) {
      if(rank == 0) {
        std::cout << "dt = " << dt << std::endl;
        time0 = comp_time_of_previous_simuls + MPI_Wtime();
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

      // Compute displacement current with Nanz method
      // Build vector
      Bsol1.get_owned_data().assign(Bsol1.local_size(),0.0);
      Bsol2.get_owned_data().assign(Bsol2.local_size(),0.0);
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
      
      // Build Matrix
      B.reset();
      bim3a_reaction (tmsh, delta0, zeta0, B, ord_displ_curr[0], ord_displ_curr[0]);
      bim3a_advection_diffusion (tmsh, sigmaB, zero_q1, B, true, ord_displ_curr[0], ord_displ_curr[1]);
      
      // Matrix * vector
      Ivec1 = B*Bsol1;
      Ivec2 = B*Bsol2;

      // Sum elements in the vectors corresponding to border nodes
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

      // Sum with results on border nodes of other processes
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
      
      // If the error on the displacement current is small enough go on otherwise half dt and repeat
      if (err_max < tol) {
        time_in_step += dt;
        sold = sold2;
        if (rank == 0 && save_error_and_comp_time) {
          // Save error and computation times
          time1 = comp_time_of_previous_simuls + MPI_Wtime();
          error_file << std::setw(20) << std::setprecision(5) << Time + time_in_step
                     << std::setw(20) << std::setprecision(7) << err_max
                     << std::setw(20) << std::setprecision(7) << time1 - time0
                     << std::setw(20) << std::setprecision(7) << time1 - start_time
                     << std::endl;
        }
        if (rank == 0 && save_displ_current) {
          // Save displacement currents
          if (!compute_2_contacts)
            I_displ_file  << std::setw(20) << Time + time_in_step
                          << std::setw(20) << I_displ1
                          << std::endl;
          else
            I_displ_file  << std::setw(20) << Time + time_in_step
                          << std::setw(20) << I_displ1
                          << std::setw(20) << I_displ2
                          << std::endl;
        }
        if (time_in_step > DT - eps) {
          Time += DT;
          ++count;
          // Save temp solution
          if (save_temp_solution && !(count % save_every_n_steps)) {
            MPI_File_open(MPI_COMM_WORLD, (temp_solution_file_name + "_" + std::to_string(count)).c_str(),
                          MPI_MODE_CREATE | MPI_MODE_WRONLY,
                          MPI_INFO_NULL, &temp_sol);
            if (rank == 0) {
              double total_time = time1 - start_time;
              MPI_File_write_at(temp_sol, 0, &count, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
              MPI_File_write_at(temp_sol, sizeof(int), &Time, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
              MPI_File_write_at(temp_sol, sizeof(int)+sizeof(double), &total_time, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
            }
            MPI_File_write_at(temp_sol, sizeof(int)+(sold.get_range_start()+2)*sizeof(double), sold.get_owned_data().data(),
                          sold.local_size(), MPI_DOUBLE, MPI_STATUS_IGNORE);
            MPI_File_close(&temp_sol);
            MPI_Barrier(MPI_COMM_WORLD);
            if (rank == 0)
              std::cout << "saved temp solution at time " + std::to_string(Time) 
                        << " and count " << count << std::endl;
            remove(last_saved_solution.c_str());
            last_saved_solution = temp_solution_file_name + "_" + std::to_string(count);
          }
          if (save_charges) {

            // Prepare support vectors
            Jz_vec.get_owned_data().assign(Jz_vec.get_owned_data().size(), 0.);
            Jz_vec.assemble(replace_op);
            Ez_eps0_vec.get_owned_data().assign(Ez_eps0_vec.get_owned_data().size(), 0.);
            Ez_eps0_vec.assemble(replace_op);
            rho_pi_k_vec.get_owned_data().assign(rho_pi_k_vec.get_owned_data().size(), 0.);
            rho_pi_k_vec.assemble(replace_op);
            
            // Wheight through apposite library function
            bim3a_boundary_mass(tmsh, 0, 5, Jz_vec, Jz_mass);
            bim3a_boundary_mass(tmsh, 0, 5, Ez_eps0_vec, Ez_eps0_mass);
            bim3a_boundary_mass(tmsh, 0, 5, rho_pi_k_vec, free_charge_mass, ord_c[0]);
            bim3a_boundary_mass(tmsh, 0, 5, rho_pi_k_vec, p1_mass, ord_c[1]);
            bim3a_boundary_mass(tmsh, 0, 5, rho_pi_k_vec, p2_mass, ord_c[2]);
            bim3a_boundary_mass(tmsh, 0, 5, rho_pi_k_vec, p3_mass, ord_c[3]);
            
            // Integrate on the border part owned by current process
            I_c = 0.;
            Ez_eps0 = 0.;
            rho_pi_k.fill(0.);
            if (Jz_vec.local_size() != rho_pi_k_vec.local_size() / (N_polcur+1) &&
                Jz_vec.local_size() != Ez_eps0_vec.local_size())
              std::cerr << "non conforming sizes" << std::endl;
            for (size_t i = 0; i < Jz_vec.local_size(); i++) {
              I_c += Jz_vec.get_owned_data()[i];
              Ez_eps0 += Ez_eps0_vec.get_owned_data()[i];
              rho_pi_k[0] += rho_pi_k_vec.get_owned_data()[i*(N_polcur+1)];
              rho_pi_k[1] += rho_pi_k_vec.get_owned_data()[i*(N_polcur+1)+1];
              rho_pi_k[2] += rho_pi_k_vec.get_owned_data()[i*(N_polcur+1)+2];
              rho_pi_k[3] += rho_pi_k_vec.get_owned_data()[i*(N_polcur+1)+3];
            }

            // Sum with that of the others
            MPI_Allreduce(MPI_IN_PLACE, &I_c, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &Ez_eps0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, rho_pi_k.data(), 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            // Print on file
            if (rank == 0) {
              if (!compute_2_contacts)
                charges_file << std::setw(20) << std::setprecision(5) << Time
                             << std::setw(20) << std::setprecision(5) << I_c
                             << std::setw(20) << std::setprecision(5) << I_displ1
                             << std::setw(20) << std::setprecision(5) << rho_pi_k[0]
                             << std::setw(20) << std::setprecision(5) << Ez_eps0 - rho_pi_k[0] + rho_pi_k[1] + rho_pi_k[2] + rho_pi_k[3]
                             << std::setw(20) << std::setprecision(5) << - rho_pi_k[1]
                             << std::setw(20) << std::setprecision(5) << - rho_pi_k[2]
                             << std::setw(20) << std::setprecision(5) << - rho_pi_k[3]
                             << std::endl;
              else
                charges_file << std::setw(20) << std::setprecision(5) << Time
                             << std::setw(20) << std::setprecision(5) << I_c
                             << std::setw(20) << std::setprecision(5) << I_displ1
                             << std::setw(20) << std::setprecision(5) << I_displ2
                             << std::setw(20) << std::setprecision(5) << rho_pi_k[0]
                             << std::setw(20) << std::setprecision(5) << Ez_eps0 - rho_pi_k[0] + rho_pi_k[1] + rho_pi_k[2] + rho_pi_k[3]
                             << std::setw(20) << std::setprecision(5) << - rho_pi_k[1]
                             << std::setw(20) << std::setprecision(5) << - rho_pi_k[2]
                             << std::setw(20) << std::setprecision(5) << - rho_pi_k[3]
                             << std::endl;
            }
          }

          // Save solution
          if (save_sol) {
            sprintf(filename, "%s%s/%s/model_1_rho_%4.4d", output_folder.c_str(), test_iter->c_str(), "sol", count);
            tmsh.octbin_export (filename, sold, ord[0]);
            sprintf(filename, "%s%s/%s/model_1_phi_%4.4d", output_folder.c_str(), test_iter->c_str(), "sol", count);
            tmsh.octbin_export (filename, sold, ord[1]);
            sprintf(filename, "%s%s/%s/model_1_p1_%4.4d",  output_folder.c_str(), test_iter->c_str(), "sol", count);
            tmsh.octbin_export (filename,sold, ord[2]);
            sprintf(filename, "%s%s/%s/model_1_p2_%4.4d",  output_folder.c_str(), test_iter->c_str(), "sol", count);
            tmsh.octbin_export (filename,sold, ord[3]);
            sprintf(filename, "%s%s/%s/model_1_p3_%4.4d",  output_folder.c_str(), test_iter->c_str(), "sol", count);
            tmsh.octbin_export (filename,sold, ord[4]);
          }

          // Enter new big time step after updating 'dt_start_big_step'
          exit_loop = true;
        }
        
        // Update dt
        dt *= std::pow(err_max/tol,-0.5)*0.9;
        dt_start_big_step = std::min(DT,std::max(dt_start_big_step*truncated_dt, dt));
        if (dt > DT - time_in_step - eps) {
          dt = DT - time_in_step;
          truncated_dt = 1;
        }
      }
      else {
        // Half dt
        dt /= 2;
        dt_start_big_step = 0.;
      }
    }
  }
  

  // Export in json format
  if (rank==0){
    error_file.close();
    charges_file.close();
    I_displ_file.close();
    if(save_charges)
      txt2json(charges_file_name, output_folder + *test_iter + "/" + "charges_file.json");
    if (save_displ_current)
      txt2json(I_displ_file_name, output_folder + *test_iter + "/" + "I_displ_file.json");
    if (save_error_and_comp_time)
      txt2json(error_file_name, output_folder + *test_iter + "/" + "error_and_comp_time.json");
  }


  // Close MPI and print report
  MPI_Barrier (MPI_COMM_WORLD);

  // Clean linear solver
  lin_solver->cleanup ();


  }
  MPI_Finalize (); 
  return 0;
}
