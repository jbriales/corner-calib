classdef CRealLidar < CScan & CConfigLidar & handle
    %CRealCamera Summary of this class goes here
    %   Detailed explanation goes here
    
    properties        
        % Old properties
        frame
%         xy % Defined by CScan
        meta
        
        % Figure handle for visualization
        hFig
        hScan
        hLines
        
        % Boolean options
        tracking
    end
    
    properties (Dependent)
        xy
    end
    
    methods
        function obj = CRealLidar( in1 )
            % obj = CRealLidar( config )
            obj = obj@CConfigLidar( in1 );
            
            obj.tracking = true;
        end
        
        function xy = get.xy( obj )
            xy = obj.v .* repmat(obj.r,2,1);
        end
        
        function setObs( obj, scan )
            obj.r  = scan.r;
            obj.ts = scan.ts;
        end
        
        function fixFrame( obj, frame )
            % Useful function to easily fix bad scan segmentation when
            % first detected in automated tracking
            obj.tracking = false;
            obj.deleteMeta;
            obj.setFrame( frame );
            obj.compute;
            obj.visualize;
            obj.tracking = true;
        end
        
        function track( obj )
            debug = 0;
            
            scan.xy = obj.xy;
            scantrack = obj.meta.scantrack;
            [scantrack, ~,~] = updateSegments( scan, scantrack, debug );
            
            obj.meta.scantrack = scantrack;
        end
        
        function [v,A_v,q,A_q] = compute( obj )
            % Takes current data loaded in Cam and tries to obtain Xi
            debug = 0;
            [v,A_v,~,~,q,A_q, lin,seg] = ...
                computeScanTO( obj.xy, obj.sd,...
                {obj.meta.scantrack.all_inliers}, debug );
            
            [obj.meta.scantrack.lin] = deal( lin{:} );
            [obj.meta.scantrack.seg] = deal( seg{:} );
            
            obj.saveMeta( obj.meta );
        end
        
        function setFrame( obj, frame )
            obj.frame = frame;
            
            if isempty(obj.frame.r)
                obj.frame.loadScan;
            end
%             obj.xy = obj.frame.xy; % Set xy points in LRF object
            obj.r = obj.frame.r; % Set xy points in LRF object
            
            % Try to load metadata
            meta = obj.loadMeta;
            if isempty(meta) % No metadata exists for current frame
                if isempty(obj.meta) ||... % The metadata property is empty for Cam
                   obj.tracking == false   % Tracking option deactivated
                    warning('No metadata available, user hint required');
                    subplot( obj.hFig )
            
                    cla
%                     delete( obj.hScan ); % Update only if necessary to improve speed
                    obj.showScan;
%                     obj.plotLIDARframe;
%                     delete( obj.hLines );
                    
%                     figure, hold on
                    meta = CMetaScan;
                    meta.manualHint( obj.xy ); % Get hint from user
                else % Case of existing previous result to do tracking
                    obj.track; % Step to update inliers
                    meta = obj.meta;
                    meta.optimized = false;
                end
                % Store new metadata in file
                obj.saveMeta( meta );
            end
            obj.meta = meta; % Set chosen metadata
        end
        function meta = loadMeta( obj )
            file = fullfile(obj.frame.metafile);
            if exist(strcat(file,'.mat'),'file')
                load( strcat(file,'.mat'), '-mat', 'meta' );
            else
%                 warning('No metadata is stored for current frame');
                meta = [];
            end
        end
        function saveMeta( obj, meta ) %#ok<INUSD>
            file = fullfile(obj.frame.metafile);
            if exist(strcat(file,'.mat'),'file')
%                 warning('Metadata is going to be overwritten');
            end
            save( strcat(file,'.mat'), 'meta' );
        end
        function deleteMeta( obj )
            file = strcat( fullfile(obj.frame.metafile), '.mat' );
            if exist(file,'file')
                delete( file );
                warning('File %s was removed',file);
            else
                warning('No metadata is stored to delete');
            end
        end
        
        function showScan( obj )
            if ~isempty(obj.frame)
                obj.hScan = plot( obj.xy(1,:), obj.xy(2,:), '.k' );
                axis equal
                hold on
                obj.plotLIDARframe
            else
                warning('No frame is loaded');
            end
        end
        
        function visualize( obj, borders )
            subplot( obj.hFig )
            cla
            
%             delete( obj.hScan ) % Update only if necessary to improve speed
            obj.showScan;
            
            if exist('borders','var')
                axis(borders);
            end
            obj.plotLIDARframe;
            
%             delete( obj.hLines )
            col = 'rgb';
            hlines = {obj.meta.scantrack.lin};
            xy = obj.xy;
            for k=1:3
                if ~isempty(hlines{k})
                    inl = cell2mat( obj.meta.scantrack(k).inliers );
%                     inl = inl{:};
                    plot( xy(1,inl), xy(2,inl), [col(k),'.'] );
                    plotHomLineWin( hlines{k}, col(k) );
                end
            end
            camroll(90); % Rotates view to put LRF pointing forwards
            
            % OLD CODE:
            % TODO: Put as option in a new class CLRF
% % %             tic
% % %             subplot(hLidar), set( gcf, 'Visible', win_visibility )
% % %             ax = axis;
% % %             cla, hold on
% % %             hold on, title('Current scan segmentation')
% % %             plot( scan.xy(1,:), scan.xy(2,:), '.k' ), axis equal
% % %             axis( ax );
% % %             plotLIDARframe( );
% % %             col = 'rgb';
% % %             for k=1:3
% % %                 if thereis_line(k)
% % %                     plot( inPts{k}(1,:), inPts{k}(2,:), [col(k),'.'] );
% % %                     plotHomLineWin( scantrack(k).lin, col(k) );
% % %                 end
% % %             end
% % %             fprintf('LIDAR PLOT TIME: %f\n',toc)
        end
        
        function plotLIDARframe( obj )
            xlabel('X [m]')
            ylabel('Y [m]')
            plot(0,0,'^', 'MarkerSize',10, 'LineWidth',3)
            plotHomLineWin( [1 0 0], 'k' );
            plotHomLineWin( [0 1 0], 'k' );
        end
        
        function h = plotPolar( obj, format )
            if ~exist('format','var')
                format = '.k';
            end
            h(1) = plot( obj.theta, obj.r, format );
        end
        
        function h = plotPts( obj, format )
            if ~exist('format','var')
                format = '.k';
            end
            h(1) = plot( obj.xy(1,:), obj.xy(2,:), format );
        end

    end
    
end

