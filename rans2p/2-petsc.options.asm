-rans2p_ksp_type                          gmres 
-rans2p_pc_type                           asm 
-rans2p_pc_asm_type                       basic 
-rans2p_ksp_max_it                        2000
-rans2p_ksp_gmres_modifiedgramschmidt 
-rans2p_ksp_gmres_restart                 300 
-rans2p_sub_ksp_type                      preonly 
-rans2p_sub_pc_factor_mat_solver_package  superlu 
-rans2p_ksp_knoll
-rans2p_sub_pc_type                       lu

-rans3p_ksp_type       gmres 
-rans3p_pc_type        hypre 
-rans3p_pc_hypre_type  boomeramg 
-rans3p_ksp_max_it     2000			

-ncls_ksp_type           gmres 
-ncls_pc_type            hypre 
-ncls_pc_hypre_type      boomeramg 
-ncls_ksp_gmres_restart  300 
-ncls_ksp_knoll 
-ncls_ksp_max_it         2000

-vof_ksp_type           gmres 
-vof_pc_type            hypre 
-vof_pc_hypre_type      boomeramg 
-vof_ksp_gmres_restart  300 
-vof_ksp_knoll 
-vof_ksp_max_it          2000

-rdls_ksp_type                          gmres 
-rdls_pc_type                           asm 
-rdls_pc_asm_type                       basic 
-rdls_ksp_gmres_modifiedgramschmidt 
-rdls_ksp_gmres_restart                 300 
-rdls_ksp_knoll 
-rdls_sub_ksp_type                      preonly 
-rdls_sub_pc_factor_mat_solver_package  superlu 
-rdls_sub_pc_type                       lu 
-rdls_ksp_max_it                        2000

-mcorr_ksp_type       cg 
-mcorr_pc_type        hypre 
-mcorr_pc_hypre_type  boomeramg 
-mcorr_ksp_max_it     2000

-clsvof_ksp_type       gmres 
-clsvof_pc_type        hypre 
-clsvof_pc_hypre_type  boomeramg 
-clsvof_ksp_max_it     2000

-mesh_ksp_type                          cg 
-mesh_pc_type                           asm 
-mesh_pc_asm_type                       basic 
-mesh_ksp_max_it                        2000
-mesh_sub_ksp_type                      preonly 
-mesh_sub_pc_factor_mat_solver_package  superlu 
-mesh_ksp_knoll 
-mesh_sub_pc_type                       lu

-phi_ksp_type       cg 
-phi_pc_type        hypre 
-phi_pc_hypre_type  boomeramg 

-pinit_ksp_type           gmres 
-pinit_pc_type            hypre 
-pinit_pc_hypre_type      boomeramg 
-pinit_ksp_gmres_restart  300 
-pinit_ksp_knoll 
-pinit_ksp_max_it         2000

-pressure_ksp_type           gmres 
-pressure_pc_type            hypre 
-pressure_pc_hypre_type      boomeramg 
-pressure_ksp_gmres_restart  300 
-pressure_ksp_knoll 
-pressure_ksp_max_it         2000
