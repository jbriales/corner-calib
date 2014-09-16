classdef CKwakComp < handle
    
    properties
        NonWeighted
        ConstWeighted
        Weighted
    end
    
    methods
        
        % Constructor
        function obj = CKwakComp (Nsim, cam_sd_N, scan_sd_N, N_co_N)
            obj.NonWeighted     = cell(cam_sd_N, scan_sd_N, N_co_N, Nsim);
            obj.ConstWeighted   = cell(cam_sd_N, scan_sd_N, N_co_N, Nsim);
            obj.Weighted        = cell(cam_sd_N, scan_sd_N, N_co_N, Nsim);
        end
        
        function optim( obj, cornerOptim, cam_sd_N_it, scan_sd_N_it, N_co_N_it, Nsim_it )
            
            [R_k_nw, t_k_nw] = cornerOptim.optimizeRt_NonWeighted;
            [R_k_cw, t_k_cw] = cornerOptim.optimizeRt_ConstWeighted;
            [R_kC_nw, t_kC_nw] = cornerOptim.optimizeRt_C_NonWeighted;
                        
            obj.NonWeighted(cam_sd_N_it, scan_sd_N_it, N_co_N_it, Nsim_it)   = {[R_k_nw  t_k_nw]};
            obj.ConstWeighted(cam_sd_N_it, scan_sd_N_it, N_co_N_it, Nsim_it) = {[R_k_cw  t_k_cw]};
            obj.Weighted(cam_sd_N_it, scan_sd_N_it, N_co_N_it, Nsim_it)      = {[R_kC_nw t_kC_nw]};
            
        end
        
    end
        
end