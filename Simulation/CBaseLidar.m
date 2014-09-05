classdef CBaseLidar < CPlane3D & CConfigLidar
    %CBaseLidar Base class for Lidar object
    % This base class stores configuration parameters and pose
    %   Constructor:
    %   Lidar = CBaseLidar( R, t, N, FOVd, sd, d_range )
    % This class inherits from CPlane3D and CConfigLidar
    %
    % See also CPlane3D, CConfigLidar.
    
    properties
        % Empty
    end
    
    methods
        % Constructor
        function obj = CBaseLidar( R, t, N, FOVd, sd, d_range )
            obj = obj@CPlane3D( R, t );
            obj = obj@CConfigLidar( N, FOVd, sd, d_range );
        end
        
        % TODO: Not working?
        function plot3( obj ) % Plot Lidar in 3D space
            h = plotframe( obj.T, 1, 'S', 'g' );
        end
    end
end