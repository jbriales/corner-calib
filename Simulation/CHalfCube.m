classdef CHalfCube < CPattern
    %CCube
    %   Constructor:
    %   pol = CCube( L, R, t )
    %       L is the size of cube side
    %
    %   p3D  - a 3x7 array with coordinates of points
    
    properties
        segments
        sd  % Standard Deviation (m) of random noise in vertexes
    end
    
    properties (SetAccess = protected) % (Read-only)
        % Empty
    end
    
    methods
        % Constructor
        function this = CHalfCube( L, R, t, sd )
            if ~exist('R','var')
                R = eye(3);
            end
            if ~exist('t','var')
                t = zeros(3,1);
            end
            if ~exist('L','var')
                L = 1;
            end
            if ~exist('sd','var')
                sd = 0;
            end
            this = this@CPattern( R, t );
            
            this.L = L;
            
            this.NF = 3;
            % Create faces: 1,2,3 correspond to planes perpendicular to
            % X,Y,Z axes respectively
            face_R{1} = -[ 0 1 0 ; 0 0 1 ; 1 0 0 ]';
            face_R{2} = -[ 0 0 1 ; 1 0 0 ; 0 1 0 ]';
            face_R{3} = -[ 1 0 0 ; 0 1 0 ; 0 0 1 ]';
            p = [ 0 0 ; L 0 ; L L ; 0 L ]';
            for i=1:this.NF
                this.face{i} = CPolygon( R * face_R{i}, t, p );
            end
            
            this.p3D = - [ 0 0 0 ;
                        L 0 0 ;
                        0 L 0 ;
                        0 0 L ;
                        0 L L ;
                        L 0 L ;
                        L L 0 ]';
            % Apply cube noise to vertexes
            this.p3D = this.p3D + sd * randn(3,7);
            % Transform vertexes to given cube pose
            this.p3D = makeinhomogeneous( this.T * makehomogeneous( this.p3D ) );
            
            % Create cell array of segments
            % First dimension: direction {x,y,z}
            % Second dimension: segment label {1,2,3}
            indexes_pairs = { [3,7],[1,2],[4,6];
                              [4,5],[1,3],[2,7];
                              [3,5],[1,4],[2,6] };
            this.segments = cell(3,3);
            pts = num2cell(this.p3D,[1,9]);
            for i=1:3
                for j=1:3
                    this.segments{i,j} = CSegment3D(pts{indexes_pairs{i,j}});
                end
            end
        end
        
        % Get array of projection of interest points in pattern
        function [uv_proj, uv_pixels] = getProjection( this, SimCamera )
            % Project extreme points of pattern (irreal results)
                [uv_proj, uv_pixels] = SimCamera.projectPattern( this );
        end
        
        % Get projection of lines
        function lines = getProjectionLines( this, SimCamera )
            [~, uv_pixels] = SimCamera.projectPattern( this );
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
        
        % Get all primitive configurations for Elqursh method
        function prim = elqursh_primitives( this )
            % For each direction, take each possible combination
            % and perm with every segment out of that direction
            prim = cell(3,18);
            for k=1:3
                pairs = nchoosek(1:3,2);
                ij = setdiff(1:3,k);
                perp_segs = this.segments(ij,:);
                perp_segs = perp_segs(:);
                temp = cell(3,6);
                for i_pairs = 1:3
                    for i_segs = 1:6
                        temp{i_pairs,i_segs} = ...
                            [ perp_segs{i_segs},...
                              this.segments{k,pairs(i_pairs,1)},...
                              this.segments{k,pairs(i_pairs,2)} ];
                    end
                end
                prim(k,:) = temp(:);
            end
            prim = prim(:);
        end
        
        % Get all primitive configurations for Trihedron method
        function prim = trihedron_primitives( this )
            % Take all possible combinations of directions of each
            % group, without order
            indexes = allcomb(1:3,1:3,1:3);
            Ncomb = size(indexes,1);
            
            prim = cell(1,Ncomb);
            for k=1:Ncomb
                C = num2cell(indexes(k,:));
                [i1,i2,i3] = C{:};
                prim{k} = [ this.segments{1,i1},...
                            this.segments{2,i2},...
                            this.segments{3,i3} ];
            end

        end
        
        % Pattern 3D representation
        function h = plot3( this ) % Plot trihedron in 3D space
            % Plot different color per direction
            color = 'rgb';
            h = zeros(3,3);
            for i=1:3
                for j=1:3
                    h(i,j) = this.segments{i,j}.plot3([color(i),'-']);
                end
            end
        end
        
        function h = plot_prim( ~, prim )
            h = zeros(1,3);
            for k=1:3
                h(k) = prim(k).plot3({'-k','LineWidth',2});
            end
        end
    end
    
end
