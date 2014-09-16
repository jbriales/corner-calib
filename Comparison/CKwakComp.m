classdef CKwakComp < handle & CBaseComp
    
    properties
        NonWeighted
        ConstWeighted
        Wasielewski
    end
    
    methods
        
        % Constructor
        function obj = CKwakComp (Nsim, cam_sd_N, scan_sd_N, N_co_N)
            obj.NonWeighted     = cell(cam_sd_N, scan_sd_N, N_co_N, Nsim);
            obj.ConstWeighted   = cell(cam_sd_N, scan_sd_N, N_co_N, Nsim);
            obj.Wasielewski     = cell(cam_sd_N, scan_sd_N, N_co_N, Nsim);
        end
        
        function optim( obj, cornerOptim )
            % Set indexes for current optimization
            cam_sd_N_it = obj.idx_cam_sd;
            scan_sd_N_it = obj.idx_scan_sd;
            N_co_N_it = obj.idx_N_co;
            Nsim_it = obj.idx_Nsim;
            
            [R_k_nw, t_k_nw] = cornerOptim.optimizeRt_NonWeighted;
            [R_k_cw, t_k_cw] = cornerOptim.optimizeRt_ConstWeighted;
            [R_kC_nw, t_kC_nw] = cornerOptim.optimizeRt_C_NonWeighted;
                        
            obj.NonWeighted(cam_sd_N_it, scan_sd_N_it, N_co_N_it, Nsim_it)   = {[R_k_nw  t_k_nw]};
            obj.ConstWeighted(cam_sd_N_it, scan_sd_N_it, N_co_N_it, Nsim_it) = {[R_k_cw  t_k_cw]};
            obj.Wasielewski(cam_sd_N_it, scan_sd_N_it, N_co_N_it, Nsim_it)   = {[R_kC_nw t_kC_nw]};
        end
        
    end
        
end