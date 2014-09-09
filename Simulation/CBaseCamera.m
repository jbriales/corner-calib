classdef CBaseCamera < CPose3D & CConfigCamera
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
        function obj = CBaseCamera( R, t, K, res, f, sd )
            obj = obj@CPose3D( R, t );
            obj = obj@CConfigCamera( K, res, f, sd );
        end
        
        function h = plot3_CameraFrustum( obj )
            % First coordinate is X (width) and second Y (height)
            corners = [ 1 1;
                        1 obj.res(2);
                        obj.res(1) obj.res(2);
                        obj.res(1) 1;
                        1 1 ]';
            frustum_corners = obj.f * hnormalise( obj.K \ makehomogeneous( corners ) );
            pts3D = makeinhomogeneous( obj.T * makehomogeneous( frustum_corners ) );
            h = zeros(1,5);
            h(5) = line( pts3D(1,:), pts3D(2,:), pts3D(3,:) );
            for i=1:4
                v = [ obj.t , pts3D(:,i) ];
                h(i) = line( v(1,:), v(2,:), v(3,:) );
            end
        end
    end
end