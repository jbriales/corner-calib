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
            % TODO
%             t = repmat(obj.t,1,obj.N);
%             dir3 = obj.R(:,1:2) * dir;
%             h = quiver3( t(1,:),t(2,:),t(3,:),...
%                          dir3(1,:), dir3(2,:), dir3(3,:) );
        end
    end
end