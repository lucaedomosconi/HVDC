function export_phi_rho (final_step, n_procs)

	export_tmesh_data ('model_0_u_%04d_%04d', {'model_0_u_%04d_%04d'}, {'phi'}, {}, {},
	                  "out_u", 0:final_step, n_procs, false);
	                   
	export_tmesh_data ('model_0_v_%04d_%04d', {'model_0_v_%04d_%04d'}, {'rho'}, {}, {},
	                  "out_v", 0:final_step, n_procs, false);
end
