classdef CCube < CPattern
    %CCube
    %   Constructor:
    %   pol = CCube( L, R, t )
    %       L is the size of cube side
    %
    %   p3D  - a 3x7 array with coordinates of points
    
    properties (SetAccess = protected) % (Read-only)
        % Empty
    end
    
    methods
        % Constructor
        function obj = CCube( L, R, t )
            if ~exist('R','var')
                R = eye(3);
            end
            if ~exist('t','var')
                t = zeros(3,1);
            end
            if ~exist('L','var')
                L = 1;
            end
            obj = obj@CPattern( R, t );
            
            obj.L = L;
            
            obj.NF = 3;
            % Create faces: 1,2,3 correspond to planes perpendicular to
            % X,Y,Z axes respectively
            face_R{1} = -[ 0 1 0 ; 0 0 1 ; 1 0 0 ]';
            face_R{2} = -[ 0 0 1 ; 1 0 0 ; 0 1 0 ]';
            face_R{3} = -[ 1 0 0 ; 0 1 0 ; 0 0 1 ]';
            p = [ 0 0 ; L 0 ; L L ; 0 L ]';
            for i=1:obj.NF
                obj.face{i} = CPolygon( R * face_R{i}, t, p );
            end
            
            obj.p3D = - [ 0 0 0 ;
                        L 0 0 ;
                        0 L 0 ;
                        0 0 L ;
                        L L 0 ;
                        0 L L ;
                        L 0 L ]';
            obj.p3D = makeinhomogeneous( obj.T * makehomogeneous( obj.p3D ) );
        end
        
        % Get array of projection of interest points in pattern
        function [uv_proj, uv_pixels] = getProjection( obj, SimCamera )
            % Project extreme points of pattern (irreal results)
                [uv_proj, uv_pixels] = SimCamera.projectPattern( obj );
        end
        
        % Get projection of lines
        function lines = getProjectionLines( obj, SimCamera )
            [~, uv_pixels] = SimCamera.projectPattern( obj );
            x = makehomogeneous( uv_pixels );
            lines.x = { snormalize( cross( x(:,1),x(:,2) ) ),...
                        snormalize( cross( x(:,4),x(:,7) ) ),...
                        snormalize( cross( x(:,3),x(:,5) ) )};
            lines.y = { snormalize( cross( x(:,1),x(:,3) ) ),...
                        snormalize( cross( x(:,4),x(:,6) ) ),...
                        snormalize( cross( x(:,2),x(:,5) ) )};
            lines.z = { snormalize( cross( x(:,1),x(:,4) ) ),...
                        snormalize( cross( x(:,2),x(:,7) ) ),...
                        snormalize( cross( x(:,3),x(:,6) ) )};
        end
        
        % Pattern 3D representation
        function h = plot3( obj ) % Plot trihedron in 3D space
            for i=1:3
                obj.face{i}.plot3;
            end
            h = plot3( obj.p3D(1,:), obj.p3D(2,:), obj.p3D(3,:), '*' );
        end
    end
    
end