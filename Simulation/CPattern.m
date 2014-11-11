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
        %% Constructor
        function obj = CPattern( R, t )
            obj = obj@CPose3D( R, t );
        end
        
        %% Simulation methods
        % Get cell array of scans on every polygon in pattern
        function [xy, range, angles, idxs] = getScan( obj, SimLidar )
            % The order of the output is right - left
            xy = cell(1,obj.NF);
            range = cell(1,obj.NF);
            angles = cell(1,obj.NF);
            idxs = cell(1,obj.NF);
            for i=1:obj.NF
                [xy{i}, range{i}, angles{i}, idxs{i}] = ...
                    SimLidar.scanPolygon( obj.face{i} );
            end
        end
                
        % Get 1x3 cell array with correspondences (lines and points)
%         function corresp = getCorrespondence( obj, Rig )
%             [xy, range, angles, idxs] = obj.getScan( Rig.Lidar );
%             [uv_proj, uv_pixels] = obj.getProjection( Rig.Camera );
%             % Do whatever necessary here
%             corresp = cell(1,3);
%         end
            
        
        %% Plotting functions
        
        % Plot scanning points on pattern polygons
        function h = plotScan( obj, SimLidar, color )
            if ~exist('color','var')
                color = 'k';
            end
            h = zeros(obj.NF,1);
            for i=1:obj.NF
                h_ = SimLidar.plotPolygonScan( obj.face{i}, color );
                if ~isempty(h_)
                    h(i) = h_;
                end
            end
        end
        
        h = plotScanFeatures( obj, SimLidar )
        
        % Plot projected points from pattern
%         function h = plotImage( obj, SimLidar )
%             h = zeros(obj.NF,1);
%             for i=1:obj.NF
%                 h_ = SimLidar.plotPolygonScan( obj.face{i} );
%                 if ~isempty(h_)
%                     h(i) = h_;
%                 end
%             end
%         end
        
        % Plot scene with Lidar and pattern
        function plotScene( obj, varargin )
%             figure
            hold on
            obj.plotReference;
            obj.plot3;
            for i=1:nargin-1
                sensor = varargin{i};
                switch class(sensor)
                    case 'CSimLidar'
                        color = 'g';
%                         sensor.plotFrame('S',color);
                        sensor.plotFrame(NaN,color);
                        obj.plotScan( sensor, color );
                        sensor.plot3_ScanFOV( color );
%                         obj.plotScanFeatures( sensor );
                    case 'CSimCamera'
                        color = 'r';
%                         sensor.plotFrame('C',color);
                        sensor.plotFrame(NaN,color);
                        sensor.plot3_PatternProjection( obj,color );
                        sensor.plot3_CameraFrustum( color );
                end
            end
            rotate3d on, axis equal
            view(115,45)
        end
    end
    
    methods (Abstract)
        % Get array of projection of interest points in pattern
        [uv_proj, uv_pixels] = getProjection( obj, SimCamera )
        
        % 3D representation
        h = plot3( obj ) % Plot pattern in 3D space
        
        % h = plotImage( obj, SimCamera )
    end
    
end