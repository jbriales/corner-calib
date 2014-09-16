classdef CZhangComp < handle & CBaseComp
        
    properties
        Linear
    end
    
    methods
        
        % Constructor
        function obj = CZhangComp (Nsim, cam_sd_N, scan_sd_N, N_co_N)
            obj.Linear     = cell(cam_sd_N, scan_sd_N, N_co_N, Nsim);
        end
        
        function optim( obj, checkerOptim )
            % Set indexes for current optimization
            cam_sd_N_it = obj.idx_cam_sd;
            scan_sd_N_it = obj.idx_scan_sd;
            N_co_N_it = obj.idx_N_co;
            Nsim_it = obj.idx_Nsim;
            
            [R_z,t_z] = checkerOptim.optimizeRt_Zhang;            
            obj.Linear(cam_sd_N_it, scan_sd_N_it, N_co_N_it, Nsim_it)  = {[R_z t_z]};
        end
        
    end
           
end