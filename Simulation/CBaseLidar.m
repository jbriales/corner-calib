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
        
        function pts = samplePolarGrid( obj, inc_ang_deg, inc_r )
            angles = -obj.FOVd/2 : inc_ang_deg : +obj.FOVd/2;
            ranges = obj.d_min : inc_r : obj.d_max;
            
            xx = cosd(angles)' * ranges;
            yy = sind(angles)' * ranges;
            pts = [ xx(:), yy(:) ]';
        end
        function pts = sampleCartGrid( obj, inc_x )
            [X,Y] = meshgrid( -obj.d_max : inc_x : +obj.d_max );
            pts = [X(:), Y(:)]';
            rr  = sqrt(sum(pts.^2,1));
            ang = atan2( Y, X );
            ang = ang(:)';
            
            pts( :, rr < obj.d_min | ...
                    rr > obj.d_max | ...
                    abs(ang) > obj.FOVr/2 ) = [];
        end
        
        function h = plot2_LidarRays( obj )
            z = zeros(1,obj.N);
            h = quiver( z,z, obj.dir(1,:), obj.dir(2,:) );
        end

        function h = plot3_LidarRays( obj )
            
            %         function h = plot3_LidarRays( obj, idxs, range )
%             if ~exist('idxs','var')
%                 dir = obj.dir;
%             else
%                 dir = obj.dir(:,idxs);
%             end
%             if exist('range','var')
%                 if size(range,2)~=size(dir,2)
%                     error('dir and range need to have same column number');
%                 end
%                 dir = dir .* repmat(range,2,1);
%             end
            t = repmat(obj.t,1,obj.N);
            dir3 = obj.R(:,1:2) * dir;
            h = quiver3( t(1,:),t(2,:),t(3,:),...
                         dir3(1,:), dir3(2,:), dir3(3,:) );
        end
    end
end