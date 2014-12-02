classdef CCube < CPattern
    %CCube
    %   Constructor:
    %   pol = CCube( L, R, t )
    %       L is the size of cube side
    %
    %   p3D  - a 3x7 array with coordinates of points
    
    properties (SetAccess = protected) % (Read-only)
        C_p3D
        C_ind % Cell with indexes in linearised cell
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
%             face_R{1} = -[ 0 1 0 ; 0 0 1 ; 1 0 0 ]';
%             face_R{2} = -[ 0 0 1 ; 1 0 0 ; 0 1 0 ]';
%             face_R{3} = -[ 1 0 0 ; 0 1 0 ; 0 0 1 ]';
%             p = [ 0 0 ; L 0 ; L L ; 0 L ]';
%             for i=1:obj.NF
%                 obj.face{i} = CPolygon( R * face_R{i}, t, p );
%             end
            
            obj.C_p3D = cell(2,2,2);
            for i=1:2
                for j=1:2
                    for k=1:2
                        obj.C_p3D{i,j,k}= L/2 * [ (-1)^(i+1)
                                                  (-1)^(j+1)
                                                  (-1)^(k+1) ];
                    end
                end
            end
            obj.p3D = reshape( cell2mat(obj.C_p3D), 3, [] );            
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
            K = SimCamera.K;
            x = K \ x;
            
%             lines = cell(3,2,2);
            lines.x = cell(2,2);
            lines.y = cell(2,2);
            lines.z = cell(2,2);
            % X lines
            for i=1:2
                for j=1:2
                    x1 = x(:,sub2ind([2,2,2],1,i,j));
                    x2 = x(:,sub2ind([2,2,2],2,i,j));
                    lines.x{i,j} = snormalize(cross(x1,x2));
                end
            end
            % Y lines
            for i=1:2
                for j=1:2
                    x1 = x(:,sub2ind([2,2,2],j,1,i));
                    x2 = x(:,sub2ind([2,2,2],j,2,i));
                    lines.y{i,j} = snormalize(cross(x1,x2));
                end
            end
            % Z lines
            for i=1:2
                for j=1:2
                    x1 = x(:,sub2ind([2,2,2],i,j,1));
                    x2 = x(:,sub2ind([2,2,2],i,j,2));
                    lines.z{i,j} = snormalize(cross(x1,x2));
                end
            end
            
        end
        
        % Pattern 3D representation
        function h = plot3( obj ) % Plot trihedron in 3D space
%             for i=1:3
%                 obj.face{i}.plot3;
%             end
            h = plot3( obj.p3D(1,:), obj.p3D(2,:), obj.p3D(3,:), '*' );
        end
    end
    
end