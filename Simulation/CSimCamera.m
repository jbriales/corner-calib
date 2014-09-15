classdef CSimCamera < CBaseCamera
    %CSimCamera Class for simulated Camera object
    % This class stores configuration parameters and pose and ease methods
    % to simulate projection of pattern points
    %   Constructor:
    %   Camera = CSimCamera( R, t, K, res, f, sd )
    % This class inherits from CBaseCamera
    %
    % See also CBaseCamera.
    
    properties
        % Empty
    end
    
    methods
        % Constructor
        function obj = CSimCamera( R, t, K, res, f, sd )
            obj = obj@CBaseCamera( R, t, K, res, f, sd );
        end
        
        % Simulate measurements given when projecting pattern points
        function [uv_proj, uv_pixels] = projectPattern( obj, pattern )
                        
            % Transform points to camera space and project them
            p3D = obj.T \ makehomogeneous(pattern.p3D); % TODO: Prod or division?
            p3D(end,:) = [];
            
            % Normalize points to f=1 distance
            p3D = hnormalise( p3D );
            
            uv_pixels = makeinhomogeneous( obj.K * p3D );
            
            % Apply gaussian noise to pixels
            uv_pixels = uv_pixels + obj.sd * randn(2,size(p3D,2));
            
            % Update canonical projection space with noise
            uv_proj   = obj.K \ makehomogeneous( uv_pixels );
        end
        
        function [uv_proj, uv_pixels ] = projectPatternInf( obj, pattern )
            % Project center point
            % Add center point
            p_center = obj.T \ makehomogeneous(pattern.p3D(:,1));
            p_center(end,:) = [];
            % Normalize points to f=1 distance
            p_center = hnormalise( p_center );
            
            if abs(p_center(1)) > tan(obj.FOVh/2) || ...
                    abs(p_center(2)) > tan(obj.FOVv/2)
                warning('Center out of Camera FOV');
                keyboard
            end
            center_pixels = makeinhomogeneous( obj.K * p_center );
            
            % Transform trihedron planes to World frame
            trihedron_planes = cellfun( @(x)x.plane, pattern.face,...
                'UniformOutput', false );
            trihedron_planes = cell2mat( trihedron_planes );
            
            % Compute frustum normals in World frame
            % Order of frustum planes is Up-Right-Down-Left in typical
            % frame Z forwards, X right, Y down
            % The order is the same as +X rot, +Y rot, -X rot, -Y rot
            frustum_planes  = cell(1,4);
            frustum_normals = [ 0 -1 0 ; 1 0 0 ; 0 1 0 ; -1 0 0 ]';
            rot = { RotationX(+obj.FOVv/2), RotationY(+obj.FOVh/2),...
                    RotationX(-obj.FOVv/2), RotationY(-obj.FOVh/2) };
            for i=1:4
                frustum_normals(:,i) = rot{i} * frustum_normals(:,i);
                frustum_planes{i} = obj.T' \ [ frustum_normals(:,i) ; 0 ];
            end
            
            % Find intersection of each pattern line (planes couple) with
            % frustum
            P_Cam = cell(1,3);
            for k=1:3 % Axis index
                idxs = setdiff(1:3,k);
                trihedron_planes_ = trihedron_planes(:,idxs);
                intersections = cellfun( @(x) makeinhomogeneous( null( [trihedron_planes_, x]' ) ), frustum_planes,...
                    'UniformOutput', false );
                intersections = cell2mat( intersections );
                ind_pos = find( intersections(k,:) > 0 );
                [~, ind_min] = min( intersections(k,ind_pos) );
                ind = ind_pos( ind_min );
                P_Cam{k} = hnormalise( makeinhomogeneous( obj.T \ makehomogeneous( intersections(:,ind) ) ) );
            end
            
            frustum_pixels = makeinhomogeneous( obj.K * cell2mat(P_Cam) );
            
            % Concatenate points
            uv_pixels = [ center_pixels, frustum_pixels ];
            
            % Apply gaussian noise to pixels
            uv_pixels = uv_pixels + obj.sd * randn(2,size(uv_pixels,2));
            
            % Update canonical projection space with noise
            uv_proj   = obj.K \ makehomogeneous( uv_pixels );
        end
        
        % Take points inside image when projected out of image
        function projectInside( obj, pts2D )
            % TODO: Check if point is contained in FOV and if not, take
            % intersection point in image border
        end
        
        % Plotting functions
        function h = plot2_PatternProjection( obj, pattern, color )
            if ~exist('color','var')
                color = 'k';
            end
            [~, uv_pixels] = projectPattern( obj, pattern );
            if ~isempty(uv_pixels)
                h = plot( uv_pixels(1,:), uv_pixels(2,:), ['*',color] );
                axis(obj.ax);
                axis ij
            else
                h = [];
            end
        end
        
        function h = plot3_PatternProjection( obj, pattern, color )
            if ~exist('color','var')
                color = 'k';
            end
%             if isa(pattern,'CTrihedron')
%                 [uv_proj, ~] = projectPatternInf( obj, pattern );
%             else
%                 [uv_proj, ~] = projectPattern( obj, pattern );
%             end
            [uv_proj, ~] = pattern.getProjection( obj );
            % Transform projection vectors with focal length f to lie on
            % plane to f distance
            uv_proj = obj.f * uv_proj;
            if ~isempty(uv_proj)
                pts3D = makeinhomogeneous( obj.T * makehomogeneous( uv_proj ) );
                h = plot3( pts3D(1,:), pts3D(2,:), pts3D(3,:), ['*',color] );
            else
                h = [];
            end
        end
    end
end