classdef CTrihedronComp < handle
    
    properties
        %Linear
        %DiagWeighted
        NonWeighted3D        
        Weighted3D
        NonWeighted2D
        Weighted2D
        Global
    end
    
    methods
        
        % Constructor
        function obj = CTrihedronComp (Nsim, cam_sd_N, scan_sd_N)
            %obj.Linear          = cell(cam_sd_N, scan_sd_N, Nsim);
            %obj.DiagWeighted    = cell(cam_sd_N, scan_sd_N, Nsim);
            obj.NonWeighted3D    = cell(cam_sd_N, scan_sd_N, Nsim);            
            obj.Weighted3D       = cell(cam_sd_N, scan_sd_N, Nsim);            
            obj.NonWeighted2D    = cell(cam_sd_N, scan_sd_N, Nsim);            
            obj.Weighted2D       = cell(cam_sd_N, scan_sd_N, Nsim);
            obj.Global           = cell(cam_sd_N, scan_sd_N, Nsim);
        end
        
        function optim( obj, triOptim, Rig, Nsim_it, cam_sd_N_it, scan_sd_N_it, WITHRANSAC )
            
            if WITHRANSAC
                triOptim.filterRotationRANSAC;
            end
            triOptim.disp_N_R_inliers;
            %R_c_s_dw = triOptim.optimizeRotation_DiagWeighted;
            R_c_s_nw = triOptim.optimizeRotation_NonWeighted;            
            R_c_s_w  = triOptim.optimizeRotation_Weighted;
            
            if WITHRANSAC
                triOptim.filterTranslationRANSAC( R_c_s_w ); % Should receive some estimated rotation
            end
            R0_for_t = R_c_s_w;
            triOptim.disp_N_t_inliers;
            triOptim.setInitialTranslation( Rig.t_c_s + 0.05*randn(3,1) );
            
            t_3D_nw = triOptim.optimizeTranslation_3D_NonWeighted( R0_for_t );
            t_3D_w  = triOptim.optimizeTranslation_3D_Weighted( R0_for_t );
            
            t_2D_nw = triOptim.optimizeTranslation_2D_NonWeighted( R0_for_t );
            t_2D_w = triOptim.optimizeTranslation_2D_Weighted( R0_for_t );

            [R_global, t_global] = triOptim.optimizeGlobal_Ort_3D( R_c_s_w, t_3D_w );  
                        
            %obj.DiagWeighted(cam_sd_N_it, scan_sd_N_it, Nsim_it) = {triOptim.optimizeRotation_DiagWeighted};
            obj.NonWeighted3D(cam_sd_N_it, scan_sd_N_it, Nsim_it)  = {[ R_c_s_nw t_3D_nw]};            
            obj.Weighted3D(cam_sd_N_it, scan_sd_N_it, Nsim_it)     = {[ R_c_s_w  t_3D_w]};            
            obj.NonWeighted2D(cam_sd_N_it, scan_sd_N_it, Nsim_it)  = {[ R_c_s_nw t_2D_nw]};            
            obj.Weighted2D(cam_sd_N_it, scan_sd_N_it, Nsim_it)     = {[ R_c_s_w  t_2D_w]};            
            obj.Global(cam_sd_N_it, scan_sd_N_it, Nsim_it)         = {[ R_global t_global]};     
            
            triOptim.resetRotationRANSAC();
            triOptim.resetTranslationRANSAC();
            
        end
        
    end
            
end
            