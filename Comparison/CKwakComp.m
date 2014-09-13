classdef CKwakComp < handle
    
    properties
        NonWeighted
        Weighted
    end
    
    methods
        
        % Constructor
        function obj = CKwakComp (Nsim, cam_sd_N, scan_sd_N)
            obj.NonWeighted     = cell(cam_sd_N, scan_sd_N, Nsim);
            obj.Weighted        = cell(cam_sd_N, scan_sd_N, Nsim);
        end
        
    end
        
end