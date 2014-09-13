classdef CZhangComp < handle
        
    properties
        Linear
    end
    
    methods
        
        % Constructor
        function obj = CZhangComp (Nsim, cam_sd_N, scan_sd_N)
            obj.Linear     = cell(cam_sd_N, scan_sd_N, Nsim);
        end
        
    end
           
end