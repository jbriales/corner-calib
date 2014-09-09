classdef CTrihedron < CPattern
    %CTrihedron 3D pattern formed by 3 orthogonal intersecting planes
    %   Constructor:
    %   pol = CTrihedron( R, t, L )
    %       L is the size of trihedron side
    %
    %   face - a 1x3 cell array formed by polygons on faces X,Y,Z
    %   p3D  - a 3x4 array with coordinates of central point, X,Y,Z extreme points, respectively
    
    properties (SetAccess = protected) % (Read-only)
        % Empty
    end
    
    methods
        % Constructor
        function obj = CTrihedron( R, t, L )
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
        function corresp = getCorrespondence( obj, Rig )
            [xy, range, angles, idxs] = obj.getScan( Rig.Lidar );
            [uv_proj, uv_pixels] = obj.getProjection( Rig.Camera );
            % Do whatever necessary here
            corresp = cell(1,3);
        end
        
        % Get
        function [N, c, A_co] = getCalibratedCornerData( obj, Camera )
            % [N, c, A_co] = getCalibratedCornerData( obj, Camera )
            % Input:
            %   img     - img structure
            %   debug   - show debugging
            %
            % Output:
            %   img     - img structure with the additional fields
            %       x   - 5xN array with: c point (in pixels), 3 angles of direction
            %       A_x - 5x5xN covariance matrices with the uncertainty of x
            
            % Control variables
            sigma = Camera.sd; % TODO: Square? sigma = sd^2?
            
            [~, img_pts] = obj.getProjection( Camera );
            
            % Define the x and A_x structures
            x   = zeros(5,1);
            A_x = zeros(5,5,1);
            
            % Assignment of auxiliar variables
            c   = img_pts(:,1);
            p   = img_pts(:,2:4) - repmat(c,1,3);
            a   = atan2( p(2,:), p(1,:) );
            
            % Estimation of x
            x(1:2) = c;
            x(3:5) = a;
            % Estimation of A_x
            A_x(1,1) = sigma;
            A_x(2,2) = sigma;
            
            A_pc = sigma * eye(4);
            J    = zeros(3,4);
            for i = 1:3
                % Set the num and den variables
                px = p(1,i);
                py = p(2,i);
                % Set the derivatives of atan2
                alpha = px.^2 + py.^2;
                d_atan_x = -py./alpha;
                d_atan_y =  px./alpha;
                % Calculate the Jacobian matrix
                J(i,1) =  d_atan_x;
                J(i,2) =  d_atan_y;
                J(i,3) = -d_atan_x;
                J(i,4) = -d_atan_y;
            end
            A_x(3:5,3:5) = J * A_pc * J';
            
            %% Convert to calibrated data
            K = Camera.K;
            % Compute homogeneous lines
            l = cell(1,3);
            for k=1:3
                n = [-sin(a(k)) cos(a(k))]';
                d = - n' * c;
                l{k} = [ n ; d ];
            end
            
            % Calibrate image
            % Central point and directions are converted from pixel values to m values
            J_an_a = zeros(3,3);
            for k=1:3
                n    = l{k}(1:2); % Normal in pixels space
                
                l{k} = K' * l{k};
                l{k} = l{k} / norm(l{k}(1:2));
                
                g = K(1:2,1:2)' * n; % Auxiliar variable: non-normalized normal in [m]
                n_n = g / norm(g); % Normal in [m] space
                
                % Compute the Jacobians of the normal directions
                J_an_nn = [-n_n(2) n_n(1)];
                J_nn_g  = [g(2)^2 -g(1)*g(2);-g(1)*g(2) g(1)^2]/norm(g)^3;
                J_g_n   = K(1:2,1:2)';
                J_n_a   = [-n(2) n(1)]';
                J_an_a(k,k) = J_an_nn * J_nn_g * J_g_n * J_n_a;
            end
            % Compute the Jacobians of the propagated point
            K_inv   = inv(K);
            c_n     = K \ makehomogeneous(c);
            J_cn_c  = [1/c_n(3) 0 -c_n(1)/c_n(3)^2;
                0 1/c_n(3) -c_n(2)/c_n(3)^2 ] * K_inv(1:3,1:2);
            c = makeinhomogeneous( c_n );
            
            %% Output data
            % Calibrated homogeneous lines
            L_P2 = cell2mat( l );
            
            % Minimal data for rotation
            N = L_P2(1:2,:);
            
            % Set uncertainty:
            J    = [J_cn_c zeros(2,3); zeros(3,2) J_an_a];
            % J    = [J_an_a zeros(3,2); zeros(2,3) J_cn_c];
            A_co = J * A_x * J';
        end
        
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