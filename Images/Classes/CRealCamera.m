classdef CRealCamera < CConfigCamera & handle
    %CRealCamera Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        frame
        img     % Current image
        img0    % Previous image, stored for tracking
        meta
        
        % Previous interesting data
        prev_ts % ts of previous frame
        
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
        
        function trackImage( obj )
            if obj.meta.optimized % If result is already optimized, skip this step
                return
            end
            
            if isempty(obj.img) || isempty(obj.img0)
                warning('img and img0 have to exist in Cam');
                return
            end
            par = struct( 'levels', 3,...
                          'iterations',5,...
                          'transform','euclidean' );
%             [ECCWarp]=iat_ecc(obj.img0,obj.img,par);
            [ECCWarp]=iat_ecc(obj.img,obj.img0,par);
            
            % Update tracking points
%             H = ECCWarp;
            R = ECCWarp(:,1:2);
            t = ECCWarp(:,3);
%             c_est = makeinhomogeneous( H * makehomogeneous(obj.meta.c) );
            c_est = R * obj.meta.c + t;
            for k=1:3
%                 q_est{k} = makeinhomogeneous( H * makehomogeneous(obj.meta.q{k}) );
                q_est{k} = R * obj.meta.q{k} + t;
                v_est{k} = snormalize( q_est{k} - c_est );
            end
            if 0
                col = 'rgb';
                for k=1:3
                    hold on
                    plot(obj.meta.q{k}(1), obj.meta.q{k}(2), ['.',col(k)]);
                    plot(q_est{k}(1), q_est{k}(2), ['.',col(k)]);
                end
                plot(c_est(1),c_est(2),'.k');
            end
            obj.meta.xi = Cxi( c_est, v_est{1}, v_est{2}, v_est{3} );
            obj.meta.c = c_est;
            obj.meta.q = q_est;
            
%             NX = size(obj.img,2);
%             NY = size(obj.img,1);
%             
%             % Compute the warped image and visualize the error
%             [wimageECC, supportECC] = iat_inverse_warping(obj.img0, ECCWarp, par.transform,1:NX,1:NY);
%             
%             %plot the warped image
%             figure;imshow(uint8(wimageECC)); title('Warped image by ECC', 'Fontsize', 14);
%             
%             % draw mosaic
%             ECCMosaic = iat_mosaic(obj.img,obj.img0,[ECCWarp; 0 0 1]);
%             figure;imshow(uint8(ECCMosaic));title('Mosaic after ECC','Fontsize',14);
        end
        
        function xi = computeXi( obj )
            % Takes current data loaded in Cam and tries to obtain Xi
            if obj.meta.optimized % If result is already optimized
                xi = obj.meta.xi;
            else % Optimizes TO xi from current metadata
                if ~isempty(obj.GOptimizer)
                    delete( obj.GOptimizer )
                end
%                 % Try to improve initial estimate with xi estimation
%                 xi0 = obj.meta.xi;
%                 if ~isempty(obj.prev_ts) && ~isempty(obj.meta.v_c) && ...
%                         all(isfinite(obj.meta.v_c)) % Rare cases with inf?
%                     inc_t = obj.frame.ts - obj.prev_ts;
%                     c_est = xi0.c.X + inc_t * obj.meta.v_c;
%                     for k=1:3
%                         rot = RotationZ( obj.meta.v_ang(k) * inc_t );
%                         rot = rot(1:2,1:2);
%                         v_est{k} = rot * xi0.v(k).X;
%                     end
%                     obj.meta.xi = Cxi( c_est, v_est{1}, v_est{2}, v_est{3} );
%                 end
                obj.GOptimizer = CGoptimizer( obj.frame.loadImg, obj.meta );
                obj.GOptimizer.expandLines = obj.linesToEnd;
                
                [xi,meta] = obj.GOptimizer.compute; % Compute new xi from initial estimate
                % Store xi velocity information for next frame
%                 if ~isempty(obj.prev_ts)
%                     inc_t = obj.frame.ts - obj.prev_ts;
%                     [meta.v_c, meta.v_ang] = xi.velocity( xi0, inc_t );
%                 end
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
            % Store previous frame interesting information
            if isempty(obj.frame)
                obj.prev_ts = [];
            else
                obj.prev_ts = obj.frame.ts;
            end
            
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
%                 warning('Metadata is going to be overwritten');
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
            if ~isempty(obj.img)
                obj.hImg = imshow( obj.img );
            else
                warning('No image was loaded yet');
            end
        end
        
        function loadImage( obj )
            if ~isempty(obj.frame)
                if ~isempty(obj.img)
                    obj.img0 = obj.img; % Backup previous image
                end
                obj.img = obj.frame.loadImg;
                if ~isempty(obj.distortion)
                    % If exists distortion, undistort
                    obj.img = cv.undistort( obj.img, obj.K, obj.distortion );
                end
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
        
        function SimCam = cloneAsSimCam( obj, cam_sd )
            
            if ~exist('cam_sd','var')
                cam_sd = obj.sd;
            end
            SimCam = CSimCamera( eye(3), zeros(3,1), obj.K,...
                    obj.res, obj.f, cam_sd );
        end
        
    end
end

