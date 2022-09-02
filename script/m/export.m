function export (final_step, n_procs)

	addpath(genpath('../../../script/m'))

	export_tmesh_data ('model_0_u_%04d_%04d', {'model_0_u_%04d_%04d'}, {'phi'}, {}, {},
	                  "out_u", 0:final_step, n_procs, false);
	                   
	export_tmesh_data ('model_0_v_%04d_%04d', {'model_0_v_%04d_%04d'}, {'rho'}, {}, {},
	                  "out_v", 0:final_step, n_procs, false);
end
                  
#export_tmesh_data ('cahn_hilliard_u_%04d_%04d', {}, {}, {'cahn_hilliard_u_%04d_%04d'}, {'u'},
#                   "out", 0:150, 1);
