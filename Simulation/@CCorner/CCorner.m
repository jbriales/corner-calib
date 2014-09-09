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
        function obj = CCorner( R, t, L, betad )
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
        
        % Get 2x3 cell array with correspondences (lines and points)
        function corresp = getCorrespondence( obj, Rig )
            
            % Create the output cell array
            corresp = cell(2,3);
            
            % Get the data and the parameters from the Lidar & Camera                        
            [xy, range, angles, idxs] = obj.getScan( Rig.Lidar );
            [uv_proj, uv_pixels]      = obj.getProjection( Rig.Camera );
            K       = Rig.Camera.K;
            theta   = Rig.Lidar.FOVd / Rig.Lidar.N ; 
            
            % Fit left and right scan line 
            [~, ~, ~, ~, lin_s_l] = computeSegmentSVD( xy{2} , 0, 0);
            [~, ~, ~, ~, lin_s_r] = computeSegmentSVD( xy{1} , 0, 0);
            
            % The middle point is the intersection of the two lines (Z = 0)
            M   = [lin_s_l'; lin_s_r'];
            P_c = - M(1:2,1:2) \ M(1:2,3);
            corresp{2,1} = [P_c; 0]; % Transform to 3D (Z = 0)

            % Extract the two boundary points in the scan (polar coordinates)
%             P_l_j = xy{1}(:,size(xy{1},2));
            P_l_j = xy{2}(:,size(xy{2},2));
            r     = norm(P_l_j);
            z     = atan2(P_l_j(2),P_l_j(1)) + deg2rad( theta/2 ) ;
            p_1    = [r*cos(z); r*sin(z)]; 
            lin_aux_l = [-p_1(2)/p_1(1) 1 0]';
            M   = [lin_s_l'; lin_aux_l'];
            P_l = - M(1:2,1:2) \ M(1:2,3);
            corresp{2,2} = [P_l; 0]; % Transform to 3D (Z = 0)

            P_r_j = xy{1}(:,1);
            r     = norm(P_r_j);
            z     = atan2(P_r_j(2),P_r_j(1)) - deg2rad( theta/2 ) ;
            p_2    = [r*cos(z); r*sin(z)]; 
            lin_aux_r = [-p_2(2)/p_2(1) 1 0]';
            M   = [lin_s_l'; lin_aux_r'];
            P_r = - M(1:2,1:2) \ M(1:2,3);
            corresp{2,3} = [P_r; 0]; % Transform to 3D (Z = 0)
            
            % Transform and assign the image lines to the output cell array
            for i = 1:3                
                s_p = makehomogeneous(uv_pixels(1:2,2*i-1));
                e_p = makehomogeneous(uv_pixels(1:2,2*i));
                l   = cross(s_p,e_p);   l = l./sqrt(l(1)^2+l(2)^2);
                corresp{1,i} = l;         
            end          
            
        end
        
        % Call to the optim method in an extern file
        obj = optim(obj, corresp, x0, weighted, Rig);
        
        % Pattern 3D representation
        function h = plot3( obj ) % Plot trihedron in 3D space
            for i=1:2
                obj.face{i}.plot3;
            end
            h = plot3( obj.p3D(1,:), obj.p3D(2,:), obj.p3D(3,:), '*' );
        end
    end
end