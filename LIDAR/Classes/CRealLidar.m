classdef CRealLidar < CScan & CConfigLidar & handle
    %CRealCamera Summary of this class goes here
    %   Detailed explanation goes here
    
    properties        
        % Old properties
        frame
%         xy % Defined by CScan
        meta
        segs
        
        % Segmentation parameters
        thres_d
        
        % Figure handle for visualization
        hFig
        hScan
        hLines
        
        % Boolean options
        tracking
        AUTOMATED
    end
    
    properties (Dependent)
        xy
    end
    
    methods
        function obj = CRealLidar( in1 )
            % obj = CRealLidar( config )
            obj = obj@CConfigLidar( in1 );
            
            obj.thres_d  = obj.sd; % Use sd as distance threshold
            
            obj.tracking = true;
            obj.AUTOMATED = true;
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
        
        function seg = computeSegment( obj, seg )
            % seg = computeSegment( obj, seg )
            % Computes the line parameters (included covariance) and store
            % in a segment struct
            [l, p, v, ~] = adjustLine( seg.pts );
            [A_l, A_p, A_ang] = propagateLineCov( seg.pts, obj.sd );
            seg.v = v;
            seg.p = p;
            seg.l = l;
            seg.A_v = A_ang;
            seg.A_p = A_p;
            seg.A_l = A_l;
        end
        
        function setFrame( obj, frame )
            obj.frame = frame;
            
            if isempty(obj.frame.r)
                obj.frame.loadScan;
            end
%             obj.xy = obj.frame.xy; % Set xy points in LRF object
            obj.r = obj.frame.r; % Set xy points in LRF object
            
            % Complete RANSAC
            if obj.AUTOMATED % Full-automated
                % Try to load metadata
                segs = obj.loadMeta;
                if isempty(segs)
                    pts = obj.xy;
                    label_0 = 1;
                    if isempty(obj.segs)
                        % Nothing
                    else % Case of existing previous result to do tracking
                        
                        % TODO: take points corresponding to last outliers
                        % and find close points to new directions to have a
                        % first good estimation
                        segs = obj.segs;
                        mask_empty = false(1,length(segs));
                        for i=1:length(segs)
                            % Find points corresponding to previous indexes
                            seg_pts = obj.xy( :, segs(i).inl );
                            % Adjust line and find outliers
                            l = adjustLine( seg_pts );
                            d = (l' * makehomogeneous( seg_pts ))';
                            seg_pts( :, abs(d) > obj.thres_d ) = [];
                            % Second refinement to recover most inliers
                            l = adjustLine( seg_pts );
                            d = (l' * makehomogeneous( obj.xy ))';
                            seg_pts = obj.xy(:, abs(d) < 3*obj.thres_d); % Increase threshold to make more flexible
                            % Adjust new line (with inliers only)
                            if isempty( seg_pts ) || size(seg_pts,2) < 2
                                mask_empty(i) = true;
                            else
                                segs(i).pts = seg_pts;
                                segs(i) = obj.computeSegment( segs(i) );
                            end
                        end
                        segs( mask_empty ) = [];
                        % Remove all just used points from array
                        used_inliers = [segs.inl];
                        pts(:,used_inliers) = []; % Remove used inliers
                        
                        % Set first new label for undetected segments
                        label_0 = max([segs.lab]) + 1;
                        
                        % Plotting
                        if 0
                            obj.segs = segs;
                            obj.visualize_segmentation;
                        end
%                         figure, obj.showScan; hold on;
%                         plot(segs_pts(1,:),segs_pts(2,:),'.c')
%                         plotHomLineWin( l, 'r' )
%                         figure, plot(d,'.k')
                    end
                    % Apply RANSAC on the rest of points
                    if 1
                        segs_INCR   = SegmentLinesIncremental( pts, label_0,...
                            obj.thres_d, 30, obj.sd );
                        if ~isempty(segs_INCR)
                            thereAreNewSegments = true;
                            segs = [segs, segs_INCR];
                        else
                            thereAreNewSegments = false;
                        end
                    else
                        segs_RANSAC = SegmentLinesRANSAC( pts, label_0,...
                            obj.thres_d, 30, obj.sd );
                        segs = [segs, segs_RANSAC];
                    end
                    
                    % Plotting
                    obj.segs = segs;
                    obj.visualize_segmentation;
                    % Filter found segments with median of dir covariance
                    if 0
                        m = median([segs.A_v]);
                        idx = [segs.A_v] > m; % Indexes to drop out
                        segs(idx) = [];
                    end
                    
                    % Plotting
                    obj.segs = segs;
                    obj.visualize_segmentation;
                    
                    if thereAreNewSegments
                        % Fuse similar segments
                        C = combnk(1:length(segs),2)';
                        segs1 = [segs(C(1,:)).l];
                        segs2 = [segs(C(2,:)).l];
                        cr  = cross(segs1,segs2,1);
                        dist_chi = zeros(1,size(C,2));
                        for i=1:size(C,2)
                            i1 = C(1,i);
                            i2 = C(2,i);
                            A_cr = skew(segs(i1).l) * segs(i2).A_l * skew(segs(i1).l)' + ...
                                skew(segs(i2).l) * segs(i1).A_l * skew(segs(i2).l)';
                            % Normalization to projective space of covariance
                            % ¡! Normalization is not applied here because the
                            % output entity is not a homogeneous vector but a metric
                            %                         Om = null( cr(:,i)' );
                            %                         A_cr = Om*Om' * A_cr * (Om*Om')';
                            dist_chi(i) = cr(:,i)' * pinv(A_cr) * cr(:,i);
                        end
                        %                     find( dist_chi < chi2inv(0.9999,3) ); % 95% probability of significance
                        idx_fuse = find( dist_chi < 100 ); % Empirical, refine this propabilistic framework
                        
                        if ~isempty(idx_fuse)
                            for i=1:size(idx_fuse)
                                i1 = C(1,idx_fuse(i));
                                i2 = C(2,idx_fuse(i));
                                %                         seg_fuse(i).inl = [segs(i1).inl, segs(i2).inl]; %
                                %                         Not set yet
                                seg_fuse(i) = segs(i1);
                                seg_fuse(i).pts = [segs(i1).pts, segs(i2).pts];
                                [l, p, v, n] = adjustLine( seg_fuse(i).pts );
                                seg_fuse(i).v = v;
                                seg_fuse(i).p = p;
                                seg_fuse(i).l = l;
                            end
                            % Substitute all new segments
                            segs(C(1,idx_fuse)) = seg_fuse;
                            segs(C(2,idx_fuse)) = [];
                        end
                    end
                    
                    % Reassign points to correct segment
                    segs = obj.refineSegments( segs );
                    
                    % TODO: Recompute line parameters with final set of
                    % points
                    for i=1:length(segs)
                        [l, p, v, n] = adjustLine( segs(i).pts );
                        [A_l, A_p, A_ang] = propagateLineCov( segs(i).pts, obj.sd );
                        segs(i).v = v;
                        segs(i).p = p;
                        segs(i).l = l;
                        segs(i).A_v = A_ang;
                        segs(i).A_p = A_p;
                        segs(i).A_l = A_l;
                    end
                    
                    % Plotting
                    obj.segs = segs;
                    obj.visualize_segmentation;
                    
                    % Store new metadata in file
                    obj.saveMeta( segs );
                end
                obj.segs = segs; % Set chosen metadata
            else
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

                        % figure, hold on
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
        end
        
        function segs = refineSegments( obj, segs )
            % Function that takes a set of segments (lines) and find more
            % likely points belonging to them in scan
            pts  = makehomogeneous(obj.xy);
            Nsegs = length(segs);
            D = zeros(Nsegs,size(pts,2));
            for i=1:Nsegs
                D(i,:) = segs(i).l' * pts;
            end

            [d,idx] = min( abs(D),[],1 );
            
            % Apply threshold (same as RANSAC?)
            idx( d > obj.thres_d ) = NaN;
            pts = obj.xy;
            mask_empty = false(1,Nsegs);
            for i=1:Nsegs
                inl = find(i==idx);
                mask_empty(i) = isempty(inl) | numel(inl)<2;
                segs(i).inl = inl;
                segs(i).pts = pts(:,inl);
            end
            segs( mask_empty ) = [];
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
            % figure;
            % camroll(90);
            % axis off
            % set(gcf,'color','w');
            
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
                    hHL(k) = plotHomLineWin( hlines{k}, col(k) );
                end
            end
            % Commented because it should be done only once (not in every
            % call)
%             camroll(90); % Rotates view to put LRF pointing forwards
            
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
        
        function visualize_segmentation( obj, borders )
            subplot( obj.hFig )
            cla
            obj.showScan;
            
            segs = obj.segs; %#ok<*PROP>
            for i=1:length(segs)
                plot(segs(i).pts(1,:),segs(i).pts(2,:),...
                    '.','Color',segs(i).col)
                text(segs(i).p(1),segs(i).p(2),...
                    num2str(segs(i).lab));
            end
            
            if exist('borders','var')
                axis(borders);
            end
            obj.plotLIDARframe;
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

