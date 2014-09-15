classdef CTrihedron < CPattern
    %CTrihedron 3D pattern formed by 3 orthogonal intersecting planes
    %   Constructor:
    %   pol = CTrihedron( L, R, t )
    %       L is the size of trihedron side
    %
    %   face - a 1x3 cell array formed by polygons on faces X,Y,Z
    %   p3D  - a 3x4 array with coordinates of central point, X,Y,Z extreme points, respectively
    
    properties (SetAccess = protected) % (Read-only)
        % Empty
    end
    
    methods
        % Constructor
        function obj = CTrihedron( L, R, t )
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
            face_R{1} = [ 0 1 0 ; 0 0 1 ; 1 0 0 ]';
            face_R{2} = [ 0 0 1 ; 1 0 0 ; 0 1 0 ]';
            face_R{3} = [ 1 0 0 ; 0 1 0 ; 0 0 1 ]';
            p = [ 0 0 ; L 0 ; L L ; 0 L ]';
            for i=1:obj.NF
                obj.face{i} = CPolygon( R * face_R{i}, t, p );
            end
            
            obj.p3D = [ 0 0 0 ;
                L 0 0 ;
                0 L 0 ;
                0 0 L ]';
            obj.p3D = makeinhomogeneous( obj.T * makehomogeneous( obj.p3D ) );
        end
        
        % Get 4x3 cell array with correspondences (lines and points)
        function co = getCorrespondence( obj, Rig )
            % Camera data
            if 0 % Project points of finite pattern
                [~, img_pts] = obj.getProjection( Rig.Camera );
            else % Project intersection with lines of infinite pattern
                [~, img_pts] = Rig.Camera.projectPatternInf( obj );
            end
            [N_im, c, A_co, L_P2, A_L_P2] = obj.getCalibratedCornerData( img_pts, Rig.Camera );
            R0 = Rig.Camera.R'; % Initial estimate for R_c_w
            [R_c_w, A_R_c_w, A_eps_c_w] = obj.getWorldNormals( R0, N_im, c, A_co );
            
            % Lidar data            
            [l,A_l,p,A_p,q,A_q, lin,seg] = ...
                obj.computeScanCorner( Rig.Lidar, 0 ); % Debug = 0
            if 1
                co = CTrihedronObservation( R_c_w, A_R_c_w, L_P2, A_L_P2,...
                    l, A_l, [], [], q, A_q );
            else
                thereis_line   = cellfun(@(x)~isempty(x), l);
                thereis_corner = cellfun(@(x)~isempty(x), q);
                if all(thereis_corner) && all(thereis_line)
                    CompleteCO = true;
                else
                    CompleteCO = false;
                end
                
                %% Store data for final optimisation
                %tic
                lab = [1 2 3];
                co.complete = CompleteCO;
                co.thereis_line = thereis_line;        % Line - rotation correspondence
                co.lab_line = lab(thereis_line);
                co.thereis_corner = thereis_corner;    % Corner - translation correspondence
                co.lab_corner = lab(thereis_corner);
                
                % For rotation optimization
                co.R_c_w   = R_c_w;     % Normal vectors to world planes (in Camera SR)
                co.A_R_c_w = A_R_c_w;
                co.l       = l;          % Direction vector of LIDAR segments (in LIDAR SR)
                co.A_l     = A_l;
                co.A_lh    = A_lh;
                
                % For translation optimization
                co.L_P2    = L_P2;
                % TODO: Compute N_repr = L_P2 / norm(L_P2)
                %             co.N_repr  = N_repr;                   % Normal vectors to reprojection planes (in Camera SR)
                co.q       = q;                        % 2D (XY) Corner points (in LIDAR SR)
                co.A_q     = A_q; 
            end
        end
        
        % Get calibrated data from corner (line normals, center and
        % covariance)
        [N, c, A_co, L_P2, A_L_P2] = getCalibratedCornerData( obj, img_pts, Camera )
        
        % Get world plane normals from Calibrated Corner Data (needs
        % initialization)
        [N, A_N, A_eps] = getWorldNormals( obj, R0, N, c, A_co )

        % Pattern 3D representation
        function h = plot3( obj ) % Plot trihedron in 3D space
            for i=1:3
                obj.face{i}.plot3;
            end
            h = plot3( obj.p3D(1,:), obj.p3D(2,:), obj.p3D(3,:), '*' );
        end
        
        function h = plotImage( obj, SimCamera )
            [uv_proj, uv_pixels] = SimCamera.projectPattern( obj );
            % TODO
            h = [];
        end
    end
    
end