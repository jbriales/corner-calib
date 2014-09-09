classdef CPolygon < CPlane3D
    %CPolygon 2D quadrilateral polygon placed in 3D space
    %   Constructor:
    %   pol = CPolygon( R, t, p )
    %       p is 2x4 array with ORDERED points
    
    properties
        p   % 2x4 array with coordinates of points relative to Polygon system
    end
    
    properties (SetAccess = private, Dependent) % (Read-only), Get method
        p3D % 3x4 array with coordinates of points from World 
    end
    
    methods
        % Constructor
        function obj = CPolygon( R, t, p )
            obj = obj@CPlane3D( R, t );
            obj.p  = p;
        end
        
        % Return mask 1 for 2D points inside the polygon
        function in_mask = isInside( obj, pts2D )
            aug_p = [ obj.p , obj.p(:,1) ];
            v  = diff( aug_p, 1, 2 );
            vn = [0 -1; 1 0] * v;
            N  = size(pts2D,2);
            cross = zeros( 4, N );
            for i=1:4
                cross(i,:) = vn(:,i)' * (pts2D - repmat(obj.p(:,i),1,N));
            end
%             cross = vn' * pts2D;
            in_mask = all( cross >= 0, 1 );
        end
                
        % 3D representation
        function p3D = get.p3D(obj)
            p3D = makeinhomogeneous( obj.M * makehomogeneous( obj.p ) );
        end
        function plot3( obj ) % Plot polygon in 3D space
            p3D = obj.p3D;
            p3D(:,end+1) = p3D(:,1); % Close polygon to plot
            plot3( p3D(1,:), p3D(2,:), p3D(3,:), 'o-' );
            quiver3( obj.t(1),obj.t(2),obj.t(3),...
                     obj.n(1),obj.n(2),obj.n(3) );
        end
    end    
end

