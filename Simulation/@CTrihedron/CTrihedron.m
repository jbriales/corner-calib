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
        
        % Get array of projection of interest points in pattern
        function [uv_proj, uv_pixels] = getProjection( obj, SimCam )
            if 0 % Project extreme points of pattern (irreal results)
                [uv_proj, uv_pixels] = SimCamera.projectPattern( obj );
            else % Simulate projection of pattern lines (as in reality)
                % Project center point
                % Add center point
                p_center = SimCam.T \ makehomogeneous(obj.p3D(:,1));
                p_center(end,:) = [];
                % Normalize points to f=1 distance
                p_center = hnormalise( p_center );
                
                if abs(p_center(1)) > tan(SimCam.FOVh/2) || ...
                        abs(p_center(2)) > tan(SimCam.FOVv/2)
                    warning('Trihedron: Center out of Camera FOV');
                    uv_proj = [];
                    uv_pixels = [];
                    return
                end
                center_pixels = makeinhomogeneous( SimCam.K * p_center );
                
                % Transform trihedron planes to World frame
                trihedron_planes = cellfun( @(x)x.plane, obj.face,...
                    'UniformOutput', false );
                trihedron_planes = cell2mat( trihedron_planes );
                
                % Compute frustum normals in World frame
                % Order of frustum planes is Up-Right-Down-Left in typical
                % frame Z forwards, X right, Y down
                % The order is the same as +X rot, +Y rot, -X rot, -Y rot
                frustum_planes  = cell(1,4);
                frustum_normals = [ 0 -1 0 ; 1 0 0 ; 0 1 0 ; -1 0 0 ]';
                rot = { RotationX(+SimCam.FOVv/2), RotationY(+SimCam.FOVh/2),...
                    RotationX(-SimCam.FOVv/2), RotationY(-SimCam.FOVh/2) };
                for i=1:4
                    frustum_normals(:,i) = rot{i} * frustum_normals(:,i);
                    frustum_planes{i} = SimCam.T' \ [ frustum_normals(:,i) ; 0 ];
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
                    P_Cam{k} = hnormalise( makeinhomogeneous( SimCam.T \ makehomogeneous( intersections(:,ind) ) ) );
                    
                    if isempty(P_Cam{k})
                        % Could happen because inf_pt is projected INSIDE FOV
                        warning('Trihedron: Not all points were projected')
                        uv_proj = [];
                        uv_pixels = [];
                        return
                    end
                end
                frustum_pixels = makeinhomogeneous( SimCam.K * cell2mat(P_Cam) );
                
                % Concatenate points
                uv_pixels = [ center_pixels, frustum_pixels ];
                
                % Apply gaussian noise to pixels
                uv_pixels = uv_pixels + SimCam.sd * randn(2,size(uv_pixels,2));
                
                % Update canonical projection space with noise
                uv_proj   = SimCam.K \ makehomogeneous( uv_pixels );
            end
        end
        
        function obj_xi = simulateTOimage( obj, img_pts, sd )
            % Control variables
            sigma = sd;
            A_cp  = sigma * eye(8);
            
            % Define the x and A_x structures
            x   = zeros(5,1);
            A_x = zeros(5,5,1);
            
            % Assignment of auxiliar variables
            c     = img_pts(:,1);
            delta = img_pts(:,2:4) - repmat(c,1,3);
            
            % Normalization of p directions
            J_v__ = cell(1,3);
            v = Manifold.S1.empty(0,3);
            for i=1:3
                v(i) = Manifold.S1( delta(:,i) );
                J_normalize = 1/norm(delta(:,i))^3 * ...
                    (eye(2) - delta(:,i)*delta(:,i)');
                J_v__{i} = v(i).DLie * J_normalize;
            end
            J = [ eye(2)            , zeros(2,3*2)      ;
                  -cell2mat(J_v__)' , blkdiag(J_v__{:}) ];
            A_xi = J * A_cp * J';
            
            % Store results in object
            obj_xi = Cxi( c, v(1), v(2), v(3) );
            obj_xi.setRepresentationCov( A_xi );
        end
        
        % Get 4x3 cell array with correspondences (lines and points)
        function co = getCorrespondence( obj, Rig )
            % Camera data
            [~, img_pts] = obj.getProjection( Rig.Camera );
            if ~isempty(img_pts)
                obj_xi  = obj.simulateTOimage( img_pts, Rig.Camera.sd );
                [obj_Nbp, obj_LP2] = obj.getBackprojectedNormals( obj_xi, Rig.Camera.K );
                R0 = Rig.Camera.R'; % Initial estimate for R_c_w
                obj_Rtri = obj.getTrihedronNormals( obj_Nbp, R0 );
            else
                co = [];
                return
            end
            
            % Lidar data            
            [l,A_l,~,~,q,A_q, ~,~] = ...
                obj.computeScanCorner( Rig.Lidar, 0 ); % Debug = 0

            co = CTrihedronObservation( obj_Rtri, obj_LP2, obj_Nbp,...
                l, A_l, [], [], q, A_q );
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