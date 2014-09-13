classdef CTrihedronComp < handle
    
    properties
        Linear
        NonWeighted
        DiagWeighted
        Weighted
    end
    
    methods
        
        % Constructor
        function obj = CTrihedronComp (Nsim, cam_sd_N, scan_sd_N)
            obj.Linear          = cell(cam_sd_N, scan_sd_N, Nsim);
            obj.NonWeighted     = cell(cam_sd_N, scan_sd_N, Nsim);
            obj.DiagWeighted    = cell(cam_sd_N, scan_sd_N, Nsim);
            obj.Weighted    = cell(cam_sd_N, scan_sd_N, Nsim);
        end
        
    end
            
end
            