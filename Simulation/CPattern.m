classdef (Abstract) CPattern < CPose3D
    %CPattern Abstract 3D pattern formed by N intersecting polygons and
    %with NP interest points
    %   Constructor is object-dependent
    % See also: CTrihedron, CCheckerboard, CCorner
    
    properties (SetAccess = protected) % (Read-only)
        L   % Pattern size
        
        NF  % Number of polygons composing the pattern
        face % 1xNF array of polygonal faces
        
        p3D % 3xNP array with coordinates of interest points from World 
    end
    
    methods
        % Constructor
        function obj = CPattern( R, t )
            obj = obj@CPose3D( R, t );
        end
        
        % Get cell array of scans on every polygon in pattern
        function [xy, range, angles, idxs] = getScan( obj, SimLidar )
            xy = cell(1,obj.NF);
            range = cell(1,obj.NF);
            angles = cell(1,obj.NF);
            idxs = cell(1,obj.NF);
            for i=1:obj.NF
                [xy{i}, range{i}, angles{i}, idxs{i}] = ...
                    SimLidar.scanPolygon( obj.face{i} );
            end
        end
        
        % Plot scanning points on pattern polygons
        function h = plotScan( obj, SimLidar )
            h = zeros(obj.NF,1);
            for i=1:obj.NF
                h_ = SimLidar.plotPolygonScan( obj.face{i} );
                if ~isempty(h_)
                    h(i) = h_;
                end
            end
        end
        
        % Plot scene with Lidar and pattern
        function plotScene( obj, SimLidar )
            figure, hold on
            SimLidar.plotReference;
            obj.plot3;
            SimLidar.plotFrame('S','g');
            obj.plotScan( SimLidar );
            rotate3d on, axis equal
        end
    end
    
    methods (Abstract)
        % 3D representation
        h = plot3( obj ) % Plot pattern in 3D space
    end
    
end