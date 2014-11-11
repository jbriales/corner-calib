classdef CSimLidar < CBaseLidar
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
        function obj = CSimLidar( R, t, N, FOVd, sd, d_range )
            obj = obj@CBaseLidar( R, t, N, FOVd, sd, d_range );
        end
        
        % Simulate measurements given when scanning polygon
        function [xy, range, angles, idxs] = scanPolygon( obj, polygon )
                        
            % Compute P2 intersection line of polygon plane with Lidar plane
            line = obj.M' * polygon.plane;
            % Compute R2 intersection points of Lidar rays with line
%             rays = obj.getSamplingLines;
%             xy = makeinhomogeneous( skew(line) * rays );
            xy = skew(line) * obj.lines;
            % Check if there is any point at infinity
            inf_idxs = abs( xy(3,:) ) < 1e-10; % Make more robust than eps
            xy = makeinhomogeneous( skew(line) * obj.lines );
            xy(:, inf_idxs) = NaN;
            
            % Transform R2 points in Lidar frame to R3 points in World frame
            pts3D = obj.transform2Dto3D( xy );
            % Transform R3 points in World frame to R2 points in polygon frame
            pts2D = polygon.transform3Dto2D( pts3D );
            % Check which points are inside polygon
            in_mask = polygon.isInside( pts2D );
            % Check negative direction intersections and set infinity
%             dir  = obj.getSamplingVectors;
%             dir_sign = dot( dir, xy, 1 );
            dir_sign = dot( obj.v, xy, 1 );
            in_mask( dir_sign < 0 ) = false;
            
            % Return measurements data inside polygon
            idxs  = find( in_mask );
            xy = xy( :, in_mask );
%             angles = obj.getSamplingAngles;
%             angles = angles( in_mask );
            angles = obj.theta( in_mask );
            range = sqrt( sum( xy.^2, 1 ) );
            % Add generated gaussian noise to Lidar measurements and update
            range = range + obj.sd * randn(1,length(range));
            xy = obj.v(:,in_mask) .* repmat(range,2,1);
        end
        
        % Plotting functions
        function h = plotPolygonScan( obj, polygon, color )
            if ~exist('color','var')
                color = 'k';
            end
            [xy, ~, ~, ~] = scanPolygon( obj, polygon );
            if ~isempty(xy)
                pts3D = obj.transform2Dto3D( xy );
                h = plot3( pts3D(1,:), pts3D(2,:), pts3D(3,:), ['.',color] );
            else
                h = [];
            end
        end
        
        function h = plot3_ScanFOV( obj, color )
            vectors2D = obj.v;
            vectors3D = obj.R(:,1:2) * vectors2D;
            
            tt = repmat( obj.t, 1, obj.N );
            h = quiver3( tt(1,:), tt(2,:), tt(3,:),...
                         vectors3D(1,:), vectors3D(2,:), vectors3D(3,:),...
                         0.25 * obj.frameSize,...
                         'Color', color );
%             h = quiver3( , 'Color', color )
%             pts3D = obj.transform2Dto3D( xy );
        end
    end
end