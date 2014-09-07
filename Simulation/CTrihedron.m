classdef CTrihedron
    %CTrihedron 3D pattern formed by 3 orthogonal intersecting planes
    %   Constructor:
    %   pol = CTrihedron( R, t, p )
    %       p is 2x4 array with ORDERED points
    
    properties (SetAccess = private) % (Read-only)
        R   % 3x3 rotation matrix of Trihedron seen from World
        t   % 3x1 translation vector of Trihedron seen from World
        T
        L   % Square pattern size
        
        face % 1x3 array of polygonal faces
        
        p3D % 3x4 array with coordinates of points from World 
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
            obj.R = R;
            obj.t = t;          
            obj.T = [ R t ; zeros(1,3) 1 ];
            
            obj.L = L;
            % Create faces
            face_R{1} = [ 0 1 0 ; 0 0 1 ; 1 0 0 ]';
            face_R{2} = [ 0 0 1 ; 1 0 0 ; 0 1 0 ]';
            face_R{3} = eye(3);
            p = [ 0 0 ; L 0 ; L L ; 0 L ]';
            for i=1:3
                obj.face{i} = CPolygon( R * face_R{i}, t, p );
            end
            
            obj.p3D = [ 0 0 0 ;
                        L 0 0 ;
                        0 L 0 ;
                        0 0 L ]';
            obj.p3D = makeinhomogeneous( obj.T * makehomogeneous( obj.p3D ) );
        end
        
        % Simulate Lidar measurements
        function [xy, range, angles, idxs] = getScan( obj, SimLidar )
            xy = cell(1,3);
            range = cell(1,3);
            angles = cell(1,3);
            idxs = cell(1,3);
            for i=1:3
                [xy{i}, range{i}, angles{i}, idxs{i}] = ...
                    SimLidar.scanPolygon( obj.face{i} );
            end
        end
        
        % 3D representation
        function plot3( obj ) % Plot trihedron in 3D space
            for i=1:3
                obj.face{i}.plot3;
            end
            plot3( obj.p3D(1,:), obj.p3D(2,:), obj.p3D(3,:), '*' );
        end
        
        function plotScan( obj, SimLidar )
            h = zeros(3,1);
            for i=1:3
                h_ = SimLidar.plotPolygonScan( obj.face{i} );
                if ~isempty(h_)
                    h(i) = h_;
                end
            end
        end
    end
    
end