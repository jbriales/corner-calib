classdef CTrihedronComp < handle
    
    properties
        %Linear
        NonWeighted
        DiagWeighted
        Weighted
    end
    
    methods
        
        % Constructor
        function obj = CTrihedronComp (Nsim, cam_sd_N, scan_sd_N)
            %obj.Linear          = cell(cam_sd_N, scan_sd_N, Nsim);
            obj.NonWeighted     = cell(cam_sd_N, scan_sd_N, Nsim);
            obj.DiagWeighted    = cell(cam_sd_N, scan_sd_N, Nsim);
            obj.Weighted    = cell(cam_sd_N, scan_sd_N, Nsim);
        end
        
        function optim( obj, triOptim, Nsim_it, cam_sd_N_it, scan_sd_N_it )
            
            triOptim.filterRotationRANSAC;
            %triOptim.disp_N_R_inliers;
            
            obj.NonWeighted(cam_sd_N_it, scan_sd_N_it, Nsim_it)  = {triOptim.optimizeRotation_NonWeighted};
            obj.DiagWeighted(cam_sd_N_it, scan_sd_N_it, Nsim_it) = {triOptim.optimizeRotation_DiagWeighted};
            obj.Weighted(cam_sd_N_it, scan_sd_N_it, Nsim_it)     = {triOptim.optimizeRotation_Weighted};        
            
        end
        
    end
            
end
            