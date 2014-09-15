classdef CZhangComp < handle
        
    properties
        Linear
    end
    
    methods
        
        % Constructor
        function obj = CZhangComp (Nsim, cam_sd_N, scan_sd_N)
            obj.Linear     = cell(cam_sd_N, scan_sd_N, Nsim);
        end
        
        function optim( obj, T_planes, lidar_points, Nsim_it, cam_sd_N_it, scan_sd_N_it )
          
            [T, ~,~,~,~] = lccZhang(T_planes,lidar_points);
            x_z = pose_inverse(T); x_z(1:3,4) = x_z(1:3,4)/1000;
            
            obj.Linear(cam_sd_N_it, scan_sd_N_it, Nsim_it)  = {x_z};
           
        end
        
    end
           
end