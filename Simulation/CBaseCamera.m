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
        
        function mask_out = isInsideFrustum( obj, pts_cam )
            mask_out = zeros(1,size(pts_cam,2));
            
            good_idxs = find( pts_cam(3,:) > 0 );
            pts_cam = pts_cam( :, good_idxs ); %#ok<FNDSB>
%             pos_mask = (pts_cam(3,:) > 0);
%             pts_cam(:,~pos_mask) = [];
            if isempty( pts_cam )
                warning('No valid input points')
                return
            end
            normalised = makeinhomogeneous( pts_cam );
            mask_inside = abs(normalised(1,:)) < tan( obj.FOVh/2 ) & ...
                          abs(normalised(2,:)) < tan( obj.FOVv/2 );
            good_idxs(~mask_inside) = [];
            mask_out(good_idxs) = true;
        end
        
        function h = plot3_CameraFrustum( obj, color )
            if ~exist('color','var')
                color = 'k';
            end
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