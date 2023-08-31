function export_phi_rho (final_step, n_procs)

	export_tmesh_data ('output/model_1_phi_%04d_%04d', {'output/model_1_phi_%04d_%04d'}, {'phi'}, {}, {},
	                  "output/out_phi", 0:final_step, n_procs, false);
	                   
	export_tmesh_data ('output/model_1_rho_%04d_%04d', {'output/model_1_rho_%04d_%04d'}, {'rho'}, {}, {},
	                  "output/out_rho", 0:final_step, n_procs, false);

	export_tmesh_data ('output/model_1_p1_%04d_%04d', {'output/model_1_p1_%04d_%04d'}, {'p1'}, {}, {},
	                  "output/out_p1", 0:final_step, n_procs, false);

	export_tmesh_data ('output/model_1_p2_%04d_%04d', {'output/model_1_p2_%04d_%04d'}, {'p2'}, {}, {},
	                  "output/out_p2", 0:final_step, n_procs, false);

	export_tmesh_data ('output/model_1_p3_%04d_%04d', {'output/model_1_p3_%04d_%04d'}, {'p3'}, {}, {},
	                  "output/out_p3", 0:final_step, n_procs, false);
					  
end
