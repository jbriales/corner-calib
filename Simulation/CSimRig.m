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
                 if nargin ~= 0
                     obj.R_c_s  = R_c_s;
                     obj.t_c_s  = t_c_s;
                     R_w_c = R_w_s * R_c_s';
                     t_w_c = t_w_s - R_w_c * t_c_s;
                     
                     obj.Camera = CSimCamera( R_w_c, t_w_c,...
                         K, res, f, cam_sd );
                     obj.Lidar  = CSimLidar( R_w_s, t_w_s,...
                         N, FOVd, scan_sd, d_range );
                 end
        end
        
        % Convert points
%         function pts_cam = LRF2Cam( pts_LRF )
%             pts_cam = obj.R_c_s(:,1:2) * pts_LRF + ...
%                 repmat( obj.t_c_s, 1, size(pts_LRF,2) );
%         end
        
        function [pixels, pts_cam] = projectLRFpts( obj, pts_LRF )            
            pts_cam = obj.R_c_s(:,1:2) * pts_LRF + ...
                repmat( obj.t_c_s, 1, size(pts_LRF,2) );
            mask = obj.Camera.isInsideFrustum( pts_cam );
            pts_cam(:,~mask) = []; % Remove points not seen by camera
            pixels = makeinhomogeneous( obj.Camera.K * pts_cam );
        end
        
        function err = calibError( obj, A_R, A_t, A_p )
            % Create grid of test points for calibration error metric
            % Currently polar grid, but backprojected uniform image points
            % could be used
            pts_LRF = obj.Lidar.samplePolarGrid( 10, 1 );
            [pixels, pts_cam] = obj.projectLRFpts( pts_LRF );
            
            N = size(pixels,2);
            
            % Compute covariance in image points by propagation
            R = obj.R_c_s;
            t = obj.t_c_s;
            K = obj.Camera.K;
            CA_pixel = cell(1,size(pts_cam,2));
            for i=1:N
                pt = pts_LRF(:,i);
                JR = skew( R(:,1:2)*pt )';
                Jt = eye(3);
                Jp = R(:,1:2);
                A_pcam = JR * A_R * JR' + Jt * A_t * Jt' + Jp * A_p * Jp';
                CA_pcam{i} = A_pcam;
                trA_pcam(i) = trace( A_pcam );
                
                Jpcam   = DprojectionP2( K * pts_cam(:,i) ) * K;
                CA_pixel{i} = Jpcam * A_pcam * Jpcam';
                trA_pixel(i) = trace( Jpcam * A_pcam * Jpcam' );
            end
            
            [err,id] = max( trA_pixel );
            
            % Representation
            if 0
                figure, debugPlotPts( pixels, '.k' ), hold on
                plot(pixels(1,id),pixels(2,id), 'r*');
                for i=1:N % Plot all covariance matrices for each test point
                    plotcov2( pixels(:,i), CA_pixel{i} );
                end
                obj.Camera.setImageBorder;
                axis ij;
            end
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
