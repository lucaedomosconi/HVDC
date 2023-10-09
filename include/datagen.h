#include <string>
#include <fstream>
void print_data (std::string const & ofile = "data.json") {
  std::ofstream ostream;
  ostream.open(ofile);
  std::string object_string =
R"({
	"output_location" : "output",
	"test_to_run" : 
	[
		"test1"
	],
	"test1" : {
		"physics" : {
			"epsilon_0" : 8.8542e-12,
			"physics_plugin" : "libPhysics.so",
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
			"NUM_REFINEMENTS" : 4,
			"maxlevel" : 6,
			"T" : 2000,
			"initial_dt_for_adaptive_time_step" : 0.00001,
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
		"options" : {
			"print_solution_every_n_seconds" : 1,
			"save_sol" : true,
			"compute_charges_on_border" : true,
			"save_displ_current" : true,
			"compute_2_contacts" : true,
	  		"save_error_and_comp_time" : true
		}
	}
})";
  ostream << object_string;
  ostream.close();
  return;
}