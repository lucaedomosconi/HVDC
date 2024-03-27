function export_phi_rho_p1_6 (initial_step, final_step, n_procs)

	export_tmesh_data ('model_1_phi_%04d_%04d', {'model_1_phi_%04d_%04d'}, {'phi'}, {}, {},
	                  "out_phi", initial_step:final_step, n_procs, false);
	                   
	export_tmesh_data ('model_1_rho_%04d_%04d', {'model_1_rho_%04d_%04d'}, {'rho'}, {}, {},
	                  "out_rho", initial_step:final_step, n_procs, false);

	export_tmesh_data ('model_1_p1_%04d_%04d', {'model_1_p1_%04d_%04d'}, {'p1'}, {}, {},
	                  "out_p1", initial_step:final_step, n_procs, false);

	export_tmesh_data ('model_1_p2_%04d_%04d', {'model_1_p2_%04d_%04d'}, {'p2'}, {}, {},
	                  "out_p2", initial_step:final_step, n_procs, false);

	export_tmesh_data ('model_1_p3_%04d_%04d', {'model_1_p3_%04d_%04d'}, {'p3'}, {}, {},
	                  "out_p3", initial_step:final_step, n_procs, false);
					  
	export_tmesh_data ('model_1_p4_%04d_%04d', {'model_1_p4_%04d_%04d'}, {'p4'}, {}, {},
	                  "out_p4", initial_step:final_step, n_procs, false);

	export_tmesh_data ('model_1_p5_%04d_%04d', {'model_1_p5_%04d_%04d'}, {'p5'}, {}, {},
	                  "out_p5", initial_step:final_step, n_procs, false);

	export_tmesh_data ('model_1_p6_%04d_%04d', {'model_1_p6_%04d_%04d'}, {'p6'}, {}, {},
	                  "out_p6", initial_step:final_step, n_procs, false);
					  
end
