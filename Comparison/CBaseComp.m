classdef CBaseComp < handle
    
    properties
        idx_cam_sd  
        idx_scan_sd 
        idx_N_co    
        idx_Nsim    
    end
    
    methods
        function obj = setIndexes( obj, idx_cam_sd, idx_scan_sd, idx_N_co, idx_Nsim )
            obj.idx_cam_sd  = idx_cam_sd;
            obj.idx_scan_sd = idx_scan_sd;
            obj.idx_N_co    = idx_N_co;
            obj.idx_Nsim    = idx_Nsim;
        end
    end
    
end