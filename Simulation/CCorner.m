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
        
        % Get 1x3 cell array with correspondences (lines and points)
        function corresp = getCorrespondence( obj, Rig )
            [xy, range, angles, idxs] = obj.getScan( Rig.Lidar );
            [uv_proj, uv_pixels] = obj.getProjection( Rig.Camera );
            % Do whatever necessary here
            corresp = cell(1,3);
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