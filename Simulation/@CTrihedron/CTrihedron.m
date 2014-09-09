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
        
        % Get calibrated data from corner (line normals, center and
        % covariance)
        [N, c, A_co] = getCalibratedCornerData( obj, Camera )

        
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