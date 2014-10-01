classdef CRealCamera < CConfigCamera & handle
    %CRealCamera Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        frame
        img
        meta
        
        % Figure handle for visualization
        hFig
        hImg
        hLines
        
        % Optional sub-object: GOptimizer
        GOptimizer
        
        % Boolean options
        tracking
        linesToEnd
    end
    
    methods
        function obj = CRealCamera( in1 )
            % obj = CRealCamera( config )
            obj = obj@CConfigCamera( in1 );
            
            % Boolean variables, default values
            obj.tracking = true;
            obj.linesToEnd = [true true true];
            
            obj.GOptimizer = [];
        end
        
        function fixFrame( obj, frame )
            % Useful function to easily fix bad scan segmentation when
            % first detected in automated tracking
            obj.tracking = false;
            obj.deleteMeta;
            obj.setFrame( frame );
            obj.computeXi;
            obj.visualize;
            obj.tracking = true;
        end
        
        function xi = computeXi( obj )
            % Takes current data loaded in Cam and tries to obtain Xi
            if obj.meta.optimized % If result is already optimized
                xi = obj.meta.xi;
            else % Optimizes TO xi from current metadata
                if ~isempty(obj.GOptimizer)
                    delete( obj.GOptimizer )
                end
                obj.GOptimizer = CGoptimizer( obj.frame.loadImg, obj.meta );
                obj.GOptimizer.expandLines = obj.linesToEnd;
                
                [xi,meta] = obj.GOptimizer.compute;
                while isempty(xi) % Sth went wrong and corner is lost, manually set
                    subplot( obj.hFig )
                    delete( obj.hImg ); % Update only if necessary to improve speed
                    delete( obj.hLines );
                    obj.showImage;

                    obj.meta.manualHint;
                    obj.GOptimizer = CGoptimizer( obj.frame.loadImg, obj.meta );
                    [xi,meta] = obj.GOptimizer.compute;
                end
                obj.meta = meta;
                obj.saveMeta( meta );
            end
        end
        
        function setFrame( obj, frame )
            obj.frame = frame;
            obj.loadImage; % Find other location to save computation if not necessary
            % Try to load metadata
            meta = obj.loadMeta;
            if isempty(meta) % No metadata exists for current image
                if isempty(obj.meta) ||... % The metadata property is empty for Cam
                   obj.tracking == false   % Tracking option deactivated
                    warning('No metadata available, user hint required');
                    subplot( obj.hFig )
            
                    delete( obj.hImg ); % Update only if necessary to improve speed
                    obj.showImage;
%                     delete( obj.hLines );
                    
%                     figure, hold on
                    meta = CMetaImg;
                    meta.manualHint; % Get hint from user
                else
                    meta = obj.meta;
                    meta.optimized = false;
                end
                % Store new metadata in file
                obj.saveMeta( meta );
            end
            obj.meta = meta; % Set chosen metadata
        end
        
        function meta = loadMeta( obj )
            file = fullfile(obj.frame.path,'meta_img',obj.frame.metafile);
            if exist(file,'file')
                load( file );
            else
%                 warning('No metadata is stored for current frame');
                meta = [];
            end
        end
        function saveMeta( obj, meta ) %#ok<INUSD>
            file = fullfile(obj.frame.path,'meta_img',obj.frame.metafile);
            if exist(file,'file')
                warning('Metadata is going to be overwritten');
            end
            save( file, 'meta' );
        end
        function deleteMeta( obj )
            file = fullfile(obj.frame.path,'meta_img',obj.frame.metafile);
            if exist(file,'file')
                delete( file );
                warning('File %s was removed',file);
            else
                warning('No metadata is stored to delete');
            end
        end
        
        function showImage( obj )
            if ~isempty(obj.frame)
                obj.hImg = imshow( obj.frame.loadImg );
            else
                warning('No frame is loaded');
            end
        end
        
        function loadImage( obj )
            if ~isempty(obj.frame)
                obj.img = obj.frame.loadImg;
            else
                warning('No frame is loaded');
            end
        end
        
        function visualize( obj )
            subplot( obj.hFig )
            
            delete( obj.hImg ) % Update only if necessary to improve speed
            obj.showImage;
            
%             delete( obj.hLines )
            col = 'rgb';
            for k=1:3
                xy = [ obj.meta.c , ...
                       obj.meta.c + obj.meta.r{k} * obj.meta.xi.v(k).X ];
                % Use exact direction v
%                 xy = [obj.meta.c , obj.meta.q{k}];
                obj.hLines(k) = line( xy(1,:), xy(2,:), 'color', col(k) );
            end
        end
        
        function visualizeTrihedronReprojection( obj )
            % visualizeRotationReprojection( c, R_c_w )
            if ~exist('format','var')
                format = '--';
            end
            
            subplot( obj.hFig )
            
            delete( obj.hImg ) % Update only if necessary to improve speed
            obj.showImage;
            
%             L = mat2cell( skew(makehomogeneous(c)) * R_c_w, 3,[1 1 1]);
%             [c,L] = uncalibrateData(c,L,K);
%             
%             % imshow( img.I ); hold on
%             plot(c(1),c(2),'om');
%             
%             rgb = 'rgb';
%             for k=1:3
%                 plotHomLineWin( L{k}, [rgb(k),format] )
%             end
        end
        
    end
end

