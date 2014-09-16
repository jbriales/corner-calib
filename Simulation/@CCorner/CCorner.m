classdef CCorner < CPattern
    %CCorner 3D pattern formed by 2 intersecting planes
    %   Constructor:
    %   pol = CCorner( R, t, L, betad )
    %       L is the size of corner side
    %       betad is the angle formed by two planes in degrees
    %   
    %   face - a 1x2 cell array formed by polygons on faces X,Y
    %   p3D  - a 3x6 array with coordinates of central line, left and right
    %   faces borders, respectively
    
    properties (SetAccess = protected) % (Read-only)
        betad   % Angle formed by two planes in degrees
        betar   % Angle formed by two planes in radians
    end
    
    methods
               
        % Constructor
        function obj = CCorner( L, R, t, betad )
            if ~exist('R','var')
%                 R = eye(3);
                R = expmap( [-1 +1 0], deg2rad(-45) );
            end
            if ~exist('t','var')
                t = zeros(3,1);
            end
            if ~exist('L','var')
                L = 1;
            end
            if ~exist('betad','var')
                betad = 120;
            end
            obj = obj@CPattern( R, t );

            obj.L = L;
            obj.betad = betad;
            obj.betar = deg2rad(betad);
            obj.NF = 2;
            % Create faces: 1,2 correspond to planes perpendicular to
            % X,Y axes respectively
            ref_R = [ 0 1 0 ; 0 0 1 ; 1 0 0 ]' * RotationY(deg2rad(-45));
            face_R{1} = ref_R * RotationY(deg2rad(+betad/2));
            face_R{2} = ref_R * RotationY(deg2rad(-betad/2));
            p = [ 0 -L ; L -L ; L L ; 0 L ]';
            for i=1:obj.NF
                obj.face{i} = CPolygon( R * face_R{i}, t, p );
            end
            
%             obj.p3D = [ 0 0 0 ;
%                         0 0 L ;
%                         L 0 0 ;
%                         L 0 L ;
%                         0 L 0 ;
%                         0 L L ]';
%             obj.p3D = makeinhomogeneous( obj.T * makehomogeneous( obj.p3D ) );
            obj.p3D = zeros(3,6);
            obj.p3D(:,1:2) = obj.face{1}.p3D(:,[1 4]);
            obj.p3D(:,3:4) = obj.face{2}.p3D(:,[2 3]);
            obj.p3D(:,5:6) = obj.face{1}.p3D(:,[2 3]);            
        end
        
        % Get array of projection of interest points in pattern
        function [uv_proj, uv_pixels] = getProjection( obj, SimCam )
            % The order of the output is Center - Right - Left
            if 0 % Project extreme points of pattern (irreal results)
                [uv_proj, uv_pixels] = SimCamera.projectPattern( obj );
            else % Simulate projection of pattern lines (as in reality)
                % Project couples of points
                P_Cam = makeinhomogeneous( SimCam.T \ makehomogeneous(obj.p3D) );
                % Check all points are in front of camera
                if any( P_Cam(3,:)<=0 )
                    warning('Corner: Point behind camera')
                    uv_proj = [];
                    uv_pixels = [];
                    return
                end
                P_Img = makeinhomogeneous( SimCam.K * P_Cam );
                P_Lin = mat2cell( P_Img, 2, [2 2 2] );
                uv_pixels = cell(1,3);
                for k=1:3 % Center, Left, Right
                    mask_inside = SimCam.isInside(P_Lin{k});
                    p = { P_Lin{k}(:,1), P_Lin{k}(:,2) };
                    if ~all(mask_inside)
                        lin = cross( makehomogeneous(p{1}), makehomogeneous(p{2}) );
                        inter_pts = makeinhomogeneous( skew(lin) * SimCam.border );
                        inter_pts = inter_pts(:,SimCam.isInside(inter_pts));
                        
                        % All intersection points could lie outside
                        if isempty(inter_pts)
                            warning('Corner: There is a line out of FOV')
                            uv_proj = [];
                            uv_pixels = [];
                            return
                        end
                        
                        if size(inter_pts,2) ~= 2
                            warning('Corner: There was some mistake obtaining FOV intersections. The observation will be dismissed');
                            keyboard
                            uv_proj = [];
                            uv_pixels = [];
                            return
                        end
                        
                        for idx=1:2
                            
                            if ~mask_inside(idx)
                                [~,ind_min] = min( sum((inter_pts - repmat(p{idx},1,2)).^2,1) );
                                p{idx} = inter_pts(:,ind_min);
                            end
                        end
                    end
                    uv_pixels{k} = cell2mat( p );
                end
                uv_pixels = cell2mat( uv_pixels );
                
                % Apply gaussian noise to pixels
                uv_pixels = uv_pixels + SimCam.sd * randn(2,size(uv_pixels,2));
                
                % Update canonical projection space with noise
                uv_proj   = SimCam.K \ makehomogeneous( uv_pixels );
            end
        end
        
        % Get 2x3 cell array with correspondences (lines and points)
        function obs = getCorrespondence( obj, Rig )
            % Get the data and the parameters from the Lidar & Camera                        
            pts = obj.getScanFeatures( Rig.Lidar );
            lines = obj.getImageFeatures( Rig.Camera );
            
            % Create the output observation object
            if isempty( pts ) || isempty( lines )
                obs = [];
            else
                obs = CCornerObservation( lines, pts );
            end
            
            % DEBUG:
%             figure, hold on
%             plot(uv_pixels(1,:),uv_pixels(2,:),'*k');
%             for i=1:3
%                 p = hnormalise( Rig.Camera.K * ( Rig.R_c_s(:,1:2) * pts{i} + Rig.t_c_s ) );
%                 plot( p(1,:),p(2,:), '.c' );
%                 plotHomLineWin( lines{i}, 'y' );
%                 lines{i}' * p
%             end
%             keyboard
%             close
        end
        
        function pts = getScanFeatures(obj, SimLRF)
            % Points are stored in a 1x3 cell array in the order:
            % 1 - Central line (intersection of lines)
            % 2 - Left line (on 'Y' plane) (interpolated)
            % 3 - Right line (on 'X' plane) (interpolated)
            
            % Initialize output
            pts = cell(1,3);
            
            % Get the data and the parameters from the Lidar & Camera                        
            [xy, range, angles, idxs] = obj.getScan( SimLRF );
            theta   = SimLRF.FOVd / SimLRF.N ; 
            
            for i=1:length(xy)
                if size(xy{i},2)<=1
                    xy{i} = [];
                end
            end
            if any( cellfun(@(x)isempty(x),xy) )
                warning('SimLRF could not get points on Corner');
                pts = [];
                return
            end
            
            % Fit left and right scan line
            [~, ~, ~, ~, lin_s_l] = computeSegmentSVD( xy{2} , 0, 0);
            [~, ~, ~, ~, lin_s_r] = computeSegmentSVD( xy{1} , 0, 0);
            
            % The middle point is the intersection of the two lines (Z = 0)
            M   = [lin_s_l'; lin_s_r'];
            P_c = - M(1:2,1:2) \ M(1:2,3);
            pts{1} = P_c; % Transform to 3D (Z = 0)

            % Extract the two boundary points in the scan (polar coordinates)
            P_l_j = xy{2}(:,end);
            r     = norm(P_l_j);
            z     = atan2(P_l_j(2),P_l_j(1)) + deg2rad( theta/2 ) ;
            p_1    = [r*cos(z); r*sin(z)]; 
            lin_aux_l = [-p_1(2)/p_1(1) 1 0]';
            M   = [lin_s_l'; lin_aux_l'];
            P_l = - M(1:2,1:2) \ M(1:2,3);
            pts{2} = P_l; % Transform to 3D (Z = 0)

            P_r_j = xy{1}(:,1);
            r     = norm(P_r_j);
            z     = atan2(P_r_j(2),P_r_j(1)) - deg2rad( theta/2 ) ;
            p_2    = [r*cos(z); r*sin(z)]; 
            lin_aux_r = [-p_2(2)/p_2(1) 1 0]';
            M   = [lin_s_r'; lin_aux_r'];
            P_r = - M(1:2,1:2) \ M(1:2,3);
            pts{3} = P_r; % Transform to 3D (Z = 0)
        end
        
        function lines = getImageFeatures(obj, SimCam)
            % Lines are stored in a 1x3 cell array in the order:
            % 1 - Central line
            % 2 - Left line (on 'Y' plane)
            % 3 - Right line (on 'X' plane)
            
            % Initialize output
            lines = cell(1,3);
            
            % Get the data and the parameters from the Camera                        
            [uv_proj, uv_pixels] = obj.getProjection( SimCam );
            if isempty(uv_proj)
                lines = [];
                return
            end
            K = SimCam.K;
            % Transform and assign the image lines to the output cell array
            for i = 1:3                
                s_p = makehomogeneous(uv_pixels(1:2,2*i-1));
                e_p = makehomogeneous(uv_pixels(1:2,2*i));
                l   = cross(s_p,e_p);   l = l./sqrt(l(1)^2+l(2)^2);
                lines{i} = l;         
            end         
        end
        
        % Plot features
        function h = plotScanFeatures( obj, SimLRF )
            pts = obj.getScanFeatures( SimLRF );
            pts = cell2mat( pts );
            xy = pts(1:2,:); % Discard zeros line in LRF frame
            
            if ~exist('color','var')
                color = 'k';
            end
            if ~isempty(xy)
                pts3D = SimLRF.transform2Dto3D( xy );
                h = plot3( pts3D(1,:), pts3D(2,:), pts3D(3,:),...
                           ['o',color], 'LineWidth', 3 );
            else
                h = [];
            end
        end
        
        % Pattern 3D representation
        function h = plot3( obj ) % Plot trihedron in 3D space
            for i=1:2
                obj.face{i}.plot3;
            end
            h = plot3( obj.p3D(1,:), obj.p3D(2,:), obj.p3D(3,:), '*' );
        end
    end
end