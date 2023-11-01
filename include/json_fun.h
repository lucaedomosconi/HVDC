#ifndef JSON_FUNCTIONS_H
#define JSON_FUNCTIONS_H
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

void print_data (std::string const & ofile) {
  std::ofstream ostream;
  ostream.open(ofile);
  std::string object_string =
R"({
	"output_location" : "output",
	"test_to_run" : 
	[
		"a_test"
	],
	"a_test" : {
		"physics" : {
			"epsilon_0" : 8.8542e-12,
			"physics_plugin" : "libPhysics.so",
			"plugin_test_index" : "test1",
			"plugin_params" : {
				"tau_p1" : 1.0,
				"tau_p2" : 10.0,
				"tau_p3" : 100.0,
				"epsilon_inf_1" : 2,
				"csi1" : 0.5,
				"csi2" : 1.0,
				"csi3" : 2.0,
				"sigma" : 3.21e-14
			}
		},
		"algorithm" : {
			"num_refinements" : 4,
			"maxlevel" : 6,
			"T" : 2000,
			"initial_dt_for_adaptive_time_step" : 0.00001,
			"biggest_time_step" : 1,
			"tol_of_adaptive_time_step" : 0.005,
			"voltage_plugin" : "libVoltages.so",
			"voltage_name" : "voltage1",
			"voltage_plugin_params" : {
				"T_discharge" : 1000,
				"tau" : 2.0
			},
			"start_from_solution" : false,
			"save_temp_solution" : true,
			"temp_sol" : {
				"file_of_starting_sol" : "test1_temp_sol_240",
				"save_every_n_steps" : 10
			}
		},
		"output" : {
			"save_sol" : true,
			"compute_charges_on_border" : true,
			"save_cond_current" : true,
			"save_error_and_comp_time" : true
		}
	}
})";
  ostream << object_string;
  ostream.close();
  return;
}

void json_export (std::ifstream &is,
                  std::ofstream &os) {
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

void txt2json (std::string const & ifile,
               std::string const & ofile) {
      std::ifstream is;
      std::ofstream os;
      is.open(ifile);
      os.open(ofile);
      json_export(is, os);
      is.close();
      os.close();
    
}

int check_overwritings (int rank,
                        json & data,
                        std::string & output_folder,
                        std::vector<std::string> & test_name) {
  bool  start_from_solution,
        save_sol,
        compute_charges_on_border,
        save_cond_current,
        save_error_and_comp_time;
  int give_warning;
  std::vector<std::string> files;
  if (rank == 0) {
    for (auto test_iter = test_name.cbegin(); test_iter != test_name.cend(); ++test_iter) {
      try {start_from_solution = data[*test_iter]["algorithm"]["start_from_solution"];}
      catch (...) {std::cerr << "Error: Unable to read ["+*test_iter+"][algorithm][start_from_solution]" << std::endl; throw;}
      try {save_sol = data[*test_iter]["output"]["save_sol"];}
      catch (...) {std::cerr << "Error: Unable to read ["+*test_iter+"][output][save_sol]" << std::endl; throw;}
      try {compute_charges_on_border = data[*test_iter]["output"]["compute_charges_on_border"];}
      catch (...) {std::cerr << "Error: Unable to read ["+*test_iter+"][output][compute_charges_on_border]" << std::endl; throw;}
      try {save_cond_current = data[*test_iter]["output"]["save_cond_current"];}
      catch (...) {std::cerr << "Error: Unable to read ["+*test_iter+"][output][save_cond_current]" << std::endl; throw;}
      try {save_error_and_comp_time = data[*test_iter]["output"]["save_error_and_comp_time"];}
      catch (...) {std::cerr << "Error: Unable to read ["+*test_iter+"][output][save_error_and_comp_time]" << std::endl; throw;}
      if (!start_from_solution) {
        if (save_error_and_comp_time  && std::filesystem::exists(output_folder + *test_iter + "/error_and_comp_time.txt")) {files.push_back(output_folder + *test_iter + "/error_and_comp_time.txt");}
        if (compute_charges_on_border && std::filesystem::exists(output_folder + *test_iter + "/charges_file.txt")) {files.push_back(output_folder + *test_iter + "/charges_file.txt");}
        if (save_cond_current         && std::filesystem::exists(output_folder + *test_iter + "/currents_file.txt")) {files.push_back(output_folder + *test_iter + "/currents_file.txt");}
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
      return 1;
    }
  }
  return 0;
}

void set_params (json & data,
                 std::string test_name,
                 double & T,
                 double & epsilon_0,
                 double & DT,
                 double & dt,
                 double & tol,
                 int & save_every_n_steps,
                 bool & start_from_solution,
                 bool & save_temp_solution,
                 bool & save_sol,
                 bool & save_error_and_comp_time,
                 bool & save_charges,
                 bool & save_cond_current,
                 std::string & temp_solution_file_name
                 ) {
  try {T = data[test_name]["algorithm"]["T"];}
  catch (...) {std::cerr << "Error: Unable to read ["+test_name+"][algorithm][T]" << std::endl; throw;}
  try {start_from_solution = data[test_name]["algorithm"]["start_from_solution"];}
  catch (...) {std::cerr << "Error: Unable to read ["+test_name+"][algorithm][start_from_solution]" << std::endl; throw;}
  try {epsilon_0 = data[test_name]["physics"]["epsilon_0"];}
  catch (...) {std::cerr << "Error: Unable to read ["+test_name+"][algorithm][epsilon_0]" << std::endl; throw;}
  try {save_temp_solution = data[test_name]["algorithm"]["save_temp_solution"];}
  catch (...) {std::cerr << "Error: Unable to read ["+test_name+"][algorithm][save_temp_solution]" << std::endl; throw;}
  try {
    if (start_from_solution) {
      temp_solution_file_name = data[test_name]["algorithm"]["temp_sol"]["file_of_starting_sol"];
      if (!std::filesystem::exists(temp_solution_file_name)) {
        std::cerr << "Error: Temporary solution file \"" << temp_solution_file_name << " does not exist!" << std::endl;
        std::exit(1);
      }
    }
  }
  catch (...) {std::cerr << "Error: Unable to read ["+test_name+"][algorithm][temp_sol][file_of_starting_sol]" << std::endl; throw;}
  
  try {
    if (save_temp_solution) {
      save_every_n_steps = data[test_name]["algorithm"]["temp_sol"]["save_every_n_steps"];
    }
  }
  catch(...) {std::cerr << "Error: Unable to read ["+test_name+"][algorithm][temp_sol][save_every_n_steps]" << std::endl; throw;}
  try {dt = data[test_name]["algorithm"]["initial_dt_for_adaptive_time_step"];}
  catch (...) {std::cerr << "Error: Unable to read ["+test_name+"][algorithm][initial_dt_for_adaptive_time_step]" << std::endl; throw;}
  try {tol = data[test_name]["algorithm"]["tol_of_adaptive_time_step"];}
  catch (...) {std::cerr << "Error: Unable to read ["+test_name+"][algorithm][tol_of_adaptive_time_step]" << std::endl; throw;}
  // Set output preferences
  try {DT = data[test_name]["algorithm"]["biggest_time_step"];}
  catch (...) {std::cerr << "Error: Unable to read ["+test_name+"][algorithm][biggest_time_step]" << std::endl; throw;}
  try {save_sol = data[test_name]["output"]["save_sol"];}
  catch (...) {std::cerr << "Error: Unable to read ["+test_name+"][output][save_sol]" << std::endl; throw;}
  try {save_error_and_comp_time = data[test_name]["output"]["save_error_and_comp_time"];}
  catch (...) {std::cerr << "Error: Unable to read ["+test_name+"][output][save_error_and_comp_time]" << std::endl; throw;}
  try {save_charges = data[test_name]["output"]["compute_charges_on_border"];}
  catch (...) {std::cerr << "Error: Unable to read ["+test_name+"][output][compute_charges_on_border]" << std::endl; throw;}
  try {save_cond_current = data[test_name]["output"]["save_cond_current"];}
  catch (...) {std::cerr << "Error: Unable to read ["+test_name+"][output][save_cond_current]" << std::endl; throw;}
}

#endif