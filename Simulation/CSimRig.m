classdef CSimRig < handle
    %CSimRig Class for simulated Rig object with Camera and Lidar
    % This class stores configuration parameters and pose and ease methods
    % to simulate Lidar and Camera measurements of a pattern
    %   Constructor:
    %   Rig = CSimRig( R_w_s, t_w_s, R_c_s, t_c_s,... % Extrinsic options
    %                  N, FOVd, scan_sd, d_range,... % Lidar options
    %                  K, res, f, cam_sd ); % Camera options
    %   Important: The given absolute pose is for Lidar from World
    % See also: CSimLidar, CCamLidar.
    
    properties
        Camera  % Simulated Camera object
        Lidar   % Simulated Lidar object
        R_c_s   % Relative rotation of Lidar seen from Camera
        t_c_s   % Relative translation of Lidar seen from Camera
    end
    
    properties (Dependent)
        Rt_c_s
    end
    
    methods
        % Constructor
        function obj = ...
            CSimRig( R_w_s, t_w_s, R_c_s, t_c_s,... % Extrinsic options
                     N, FOVd, scan_sd, d_range,... % Lidar options
                     K, res, f, cam_sd ) % Camera options
            obj.R_c_s  = R_c_s;
            obj.t_c_s  = t_c_s;
            R_w_c = R_w_s * R_c_s';
            t_w_c = t_w_s - R_w_c * t_c_s;
            
            obj.Camera = CSimCamera( R_w_c, t_w_c,...
                                     K, res, f, cam_sd );
            obj.Lidar  = CSimLidar( R_w_s, t_w_s,...
                                    N, FOVd, scan_sd, d_range ); 
        end
        
        % Update poses
        function obj = updateLRFPose( obj, R_w_s, t_w_s )
            obj.Lidar.R = R_w_s;
            obj.Lidar.t = t_w_s;
            
            R_w_c = R_w_s * obj.R_c_s';
            t_w_c = t_w_s - R_w_c * obj.t_c_s;
            obj.Camera.R = R_w_c;
            obj.Camera.t = t_w_c;
        end
        
        function obj = updateCamPose( obj, R_w_c, t_w_c )
            obj.Camera.R = R_w_c;
            obj.Camera.t = t_w_c;
            
            R_w_s = R_w_c * obj.R_c_s;
            t_w_s = t_w_c + R_w_c * obj.t_c_s;
            
            obj.Lidar.R = R_w_s;
            obj.Lidar.t = t_w_s;
        end
        
        % Get-function
        function Rt_c_s = get.Rt_c_s( obj )
            Rt_c_s = [ obj.R_c_s obj.t_c_s ];
        end
        
    end
    
end
