classdef CSimCamera < CBaseCamera
    %CSimLidar Class for simulated Lidar object
    % This class stores configuration parameters and pose and ease methods
    % to simulate intersection with polygons in 3D space
    %   Constructor:
    %   Lidar = CSimLidar( R, t, N, FOVd, sd, d_range )
    % This class inherits from CBaseLidar
    %
    % Own methods:
    % [xy, range, angles, idxs] = scanPolygon( CPolygon )
    %   Get 2D points, range, angles and indexes for scan on polygon
    %
    % See also CBaseLidar.
    
    properties
        % Empty
    end
    
    methods
        % Constructor
        function obj = CSimCamera( R, t, K, res, f, sd )
            obj = obj@CBaseCamera( R, t, K, res, f, sd );
        end
        
        % Simulate measurements given when projecting pattern points
        function [uv_proj, uv_pixels] = projectPattern( obj, pattern )
                        
            % Transform points to camera space and project them
            p3D = obj.T \ makehomogeneous(pattern.p3D); % TODO: Prod or division?
            p3D(end,:) = [];
            uv_pixels = makeinhomogeneous( obj.K * p3D );
            
            % Apply gaussian noise to pixels
            uv_pixels = uv_pixels + obj.sd * randn(2,size(p3D,2));
            
            % Update canonical projection space with noise
            uv_proj   = obj.K \ makehomogeneous( uv_pixels );
        end
        
        % Plotting functions
        function h = plot2_PatternProjection( obj, pattern )
            [~, uv_pixels] = projectPattern( obj, pattern );
            if ~isempty(uv_pixels)
                % TODO: Set ij axis?
                h = plot( uv_pixels(1,:), uv_pixels(2,:), '*' );
                axis(obj.ax);
                axis ij
            else
                h = [];
            end
        end
        
        function h = plot3_PatternProjection( obj, pattern )
            [uv_proj, ~] = projectPattern( obj, pattern );
            if ~isempty(uv_proj)
                pts3D = makeinhomogeneous( obj.T * makehomogeneous( uv_proj ) );
                h = plot3( pts3D(1,:), pts3D(2,:), pts3D(3,:), '*' );
            else
                h = [];
            end
        end
    end
end