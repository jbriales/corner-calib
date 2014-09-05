classdef CSimLidar < CBaseLidar
    %CSimLidar Class for simulated Lidar object
    % This class stores configuration parameters and pose and ease methods
    % to simulate intersection with polygons in 3D space
    %   Constructor:
    %   Lidar = CSimLidar( R, t, N, FOVd, sd, d_range )
    % This class inherits from CBaseLidar
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
        
        % Measurements from scanning polygon
        function scanPolygon( obj, polygon )
            % Compute P2 intersection line of polygon plane with Lidar plane
            % TODO: Mistake here
            % figure, hold on, plotframe(eye(4),2,'W','k'), Lidar.plot3, pol.plot3
            % figure, plot(x(1,:),x(2,:),'.r'), hold on, plotHomLineWin( line, 'g' ), N = obj.N, quiver(zeros(1,N), zeros(1,N), dir(1,:), dir(2,:))
            line = makehomogeneous( obj.Ms )' * polygon.plane;
            % Compute P2 intersection points of Lidar rays with line
            dir  = obj.getSamplingVectors; % Debug
            rays = obj.getSamplingLines;
            x = makeinhomogeneous( skew(line) * rays );
            pts3D = transform2Dto3D( obj, x );
            pts2D = transform3Dto2D( polygon, pts3D );
            in_mask = inInside( polygon, pts2D );
        end
    end
end