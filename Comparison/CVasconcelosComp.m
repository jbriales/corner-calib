classdef CVasconcelosComp < handle
        
    properties
        Linear
    end
    
    methods
        
        % Constructor
        function obj = CVasconcelosComp (Nsim, cam_sd_N, scan_sd_N)
            obj.Linear     = cell(cam_sd_N, scan_sd_N, Nsim);
        end
        
    end
           
end