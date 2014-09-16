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
        
        function mask = isInside( obj, pts_img )
            pts_img = round(pts_img); % To avoid eps-order errors
            mask = pts_img(1,:) >= 1 & ...
                   pts_img(1,:) <= obj.res(1) & ...
                   pts_img(2,:) >= 1 & ...
                   pts_img(2,:) <= obj.res(2);
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