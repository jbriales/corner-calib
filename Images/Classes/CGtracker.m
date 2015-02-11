classdef CGtracker < handle
    %CGtracker Class with methods for tracking segments in image
    %   Detailed explanation goes here
    
    properties
        img
        Gmag    % Image gradient magnitude
        Gang    % Image gradient angle
        Gdir    % Image gradient direction
        imSize
        
        % Tracked elements
        segs
        maskSegs % Activate or deactivate tracking for concrete segments
        
        % Tracking parameters
        h = 5 % Half-Width of segment domain window
        g = 10 % Half-Growth of segment in its direction (for size increase)
        
        % Control parameters
        WITH_PLOT = true
        WITH_DEB  = false
        
        % Plotting
        hFigure
        hSegs   % Handle for segment lines in figure
        hTags   % Handle for segment tags in figure
        hIm     % Handle for image data in figure
    end
    
    methods
        function this = CGtracker( img, WITH_PLOT )
            if nargin~= 0
                if ndims( img ) == 3 % Color image
                    img = rgb2gray( img );
                end
                this.img = img;
                
                % Compute image gradients once
                [this.Gmag,this.Gang] = imgradient( img );
                ang = deg2rad( this.Gang(:) );
                % Gdir is stored with linear indexing to avoid 3rd dimension
                this.Gdir = [ -cos(ang) , +sin(ang) ]'; % Check sign! Image coordinates?
                this.imSize = size( img );
                
                % Default tracking parameters
                % ...
                
                if exist('WITH_PLOT','var')
                    this.WITH_PLOT = WITH_PLOT;
                end
            end
            
            if this.WITH_PLOT
                this.hFigure = figure;
            end
        end
        
        % Access methods
        function loadImage( this, img )
            % Test dimension
%             if any( size(img) ~= this.imSize )
%                 error('Bad dimension of input image wrt previous');
%             end
            
            if ndims( img ) == 3 % Color image
                img = rgb2gray( img );
            end
            this.img = img;
            
            % Compute image gradients once
            [this.Gmag,this.Gang] = imgradient( img );
            ang = deg2rad( this.Gang(:) );
            % Gdir is stored with linear indexing to avoid 3rd dimension
            this.Gdir = [ -cos(ang) , +sin(ang) ]'; % Check sign! Image coordinates?            
            
            this.imSize = size( img );
        end
        function loadSegs( this, file )
            % load(...)
            % TEMPORAL
            load( file, '-mat' );
            if exist('temp_segs','var')
                this.segs = temp_segs;
                this.maskSegs = temp_mask;
            else
                Sload = load( file, '-mat', 'tracker' );
                this.segs = Sload.tracker.segs;
                this.maskSegs = Sload.tracker.maskSegs;
            end
        end
        function loadSegs2( this, file )
            % load(...)
            Sload = load( file, '-mat', 'tracker' );
            this.segs = Sload.tracker.segs;
            this.maskSegs = Sload.tracker.maskSegs;
%             this.segs = temp_segs;
%             this.maskSegs = temp_mask;
        end
        function loadNet( this, file )
            
        end
        function toggleSegs( this, mask )
            this.maskSegs( mask ) = ~this.maskSegs( mask );
        end
        
        % Functions for initialization and user hint
        function hint( this, num )
            % hint( this )
            % Get image data for tracking by manually clicking points
            % It should be applied after showing image in some figure
            % Outputs:

            % Create figure
            h = figure('Name','User hints');
            imshow( this.img ); hold on;
            
            count = 0;
            if ~exist('num','var')
                % Initialize empty vector
                this.segs = CSegment2D;
                
                str = '';
                while ~strcmpi(str,'y')
                    this.segs(count+1) = this.hintSegment;
                    count = count + 1;
                    
                    % Plot added segment
                    this.segs(end).plot('r');
                    
                    str = input('Input [y] to finish, or Enter to add new segment: ','s');
                end
                
                this.maskSegs = true(1,count);
            else
                % Plot previous lines and tags
                tags = cellfun(@(x)num2str(x), num2cell(1:numel(this.segs)),...
                    'UniformOutput', false);
                this.segs.plot('y', tags);
                
                for i = num
                    this.segs(i) = this.hintSegment;
                    count = count + 1;
                    % Plot added segment
                    this.segs(end).plot('r');
                    this.maskSegs(i) = true;
                end
            end
            fprintf('User hints finished. %d segments added.\n',count);
            close( h );
        end
        function addNetPts( this, num )
            % addNetPts( this )
            % Add points to existing object net property
            if ~exist('num','var')
                num = [];
            end
            
            % Create figure (GUI)
%             h = figure('Name','Correct current points and segments');
            set(this.hFigure,'Name','Correct current points and segments');
%             imshow( this.img ); hold on;
            set(this.hIm,'CData',this.img);
            % Remove previous segments and tags
            if any(this.hSegs),delete(this.hSegs),end
            if any(this.hTags),delete(this.hTags),end
            
            hPts = [];
            if ~isempty( this.segs.pts )
                hPts = plot(this.segs.pts(1,:), this.segs.pts(2,:), '*y');
            end
            tags = cellfun(@(x)num2str(x), num2cell(1:this.segs.Npts),...
                        'UniformOutput', false);
            gPts = zeros(1,this.segs.Npts);
            for i=1:this.segs.Npts
                gPts(i) = text(this.segs.pts(1,i), this.segs.pts(2,i), tags{i});
            end
            set(gPts,'FontWeight','bold',...
                     'FontSize',15);
            idxs = str2num( input('Correct existing points? [indeces]: ','s') );
            % str2num is used together with 's' option to make both [a b d]
            % and a b d inputs possible
            hPts1 = zeros(1,numel(idxs));
            counterPts = 1;
            for i=idxs
                title(sprintf('Give new position for point #%d',i));
                set( gPts(i),'Color','r','FontWeight','bold' );
                % Get point to correct i-th net point
                zoom on;
                pause;
                % Clean buffer?
                input('');
                [x, y] = getpts;
                this.segs.pts(:,i) = [x;y];
                hPts1(counterPts) = plot(x,y,'*g');
                set( gPts(i),'Color','k','FontWeight','normal' );
                counterPts = counterPts + 1;
                zoom out;
            end
            
            % Get points to add to the net
            title('Add new points?');
            str = input('Add new points [y]? ','s');
            hPts2 = [];
            if strcmpi(str,'y')
                title('Give position for new points');
                [x, y] = getpts;
                newPts  = [x' ; y'];
                hPts2 = plot( newPts(1,:), newPts(2,:), 'r*' );
                
                % Put all points together
                this.segs.pts = [this.segs.pts, newPts];
            end
            Npts = this.segs.Npts;
            pts  = this.segs.pts;
            
            % Set segments which connect marked points:
            tags = this.segs.tags( this.maskSegs );
            [this.hSegs,this.hTags] = ...
                this.segs.segs(this.maskSegs).plot('y', tags);
            hSegs2 = zeros(1,numel(num));
            counter_segs = 1;
            for k = num
                title(sprintf('Give end points for segment #%d',k));
                set(this.hFigure, 'Name',['Correct segment #',num2str(k)]);
                seg = this.hintSegment;
                % Find closest points among existing
                [~,i] = min( sum((pts - repmat(seg.p1,1,Npts)).^2,1) );
                [~,j] = min( sum((pts - repmat(seg.p2,1,Npts)).^2,1) );
                
                this.segs.C(1:2,k) = [i;j];
                
                % Plot added segment
                seg.p1 = pts(:,i);
                seg.p2 = pts(:,j);
                hSegs2(counter_segs) = seg.plot('r');
                counter_segs = counter_segs + 1;
                
                % Make hinted segment unmasked
                this.maskSegs(k) = true;
            end
            
            % Finally update the segments in net
            this.segs.updateSegs;
            
            % Close hint figure
%             close(h);

            if any(gPts),  delete( gPts ),  end
            if any(hPts),  delete( hPts  ), end
            if any(hPts1), delete( hPts1 ), end
            if any(hPts2), delete( hPts2 ), end
            delete( this.hSegs );
            if any(hSegs2), delete( hSegs2 ), end
            delete( this.hTags );

            % Plot new segments
            set(this.hFigure,'Name','Tracker');
            title('');
            tags = this.segs.tags( this.maskSegs );
            [this.hSegs,this.hTags] = ...
                this.segs.segs(this.maskSegs).plot('y', tags);
        end
        function hintNet( this )
            % hintNet( this )
            % Get image data for tracking by manually clicking points
            % It should be applied after showing image in some figure
            % Outputs:

            % Create figure
            h = figure('Name','User hints');
            imshow( this.img ); hold on;
            
            % Get points (all points contained in net)
            title('Choose corner/end points');
            [x, y] = getpts;
            pts  = [x' ; y'];
            Npts = size(pts,2);
            plot( pts(1,:), pts(2,:), 'r*' );
            
            % Set segments which connect marked points
            count = 0;
            if ~exist('num','var')
                % Initialize empty vector
                this.segs = CSegment2D;
                
                str = '';
                while ~strcmpi(str,'y')
                    title(sprintf('Choose end points for segment #%d',count+1));
                    seg = this.hintSegment;
                    % Find closest points among existing
                    [~,i] = min( sum((pts - repmat(seg.p1,1,Npts)).^2,1) );
                    [~,j] = min( sum((pts - repmat(seg.p2,1,Npts)).^2,1) );
                    
                    count = count + 1;
                    C(1:2,count) = [i;j];
                    
                    % Plot added segment
                    seg.p1 = pts(:,i);
                    seg.p2 = pts(:,j);
                    seg.plot('r');
                    
                    str = input('Input [y] to finish, or Enter to add new segment: ','s');
                end
                
                this.maskSegs = true(1,count);
                
                % Create final object
                this.segs = CNetSegments( pts, C );
%             else
%                 % Plot previous lines and tags
%                 tags = cellfun(@(x)num2str(x), num2cell(1:numel(this.segs)),...
%                     'UniformOutput', false);
%                 this.segs.plot('y', tags);
%                 
%                 for i = num
%                     this.segs(i) = this.hintSegment;
%                     count = count + 1;
%                     % Plot added segment
%                     this.segs(end).plot('r');
%                     this.maskSegs(i) = true;
%                 end
            end
            fprintf('User hints finished. %d segments added.\n',count);
            close( h );
        end
        function seg = hintSegment( ~ )
            [x, y] = getpts;
            while isempty(x) || numel(x)~=2
                warning('Input 2 points');
                [x,y] = getpts;
            end
            % Create segment object from two points
            p1 = [x(1), y(1)]';
            p2 = [x(2), y(2)]';
            seg = CSegment2D( p1, p2 );
        end
        
        % High-level methods
        function IAT_update( this, newImg )
            param = struct( 'levels', 3,...
                          'iterations',5,...
                          'transform','euclidean' );
%             [ECCWarp]=iat_ecc(obj.img0,obj.img,par);
            ECCWarp = iat_ecc(newImg, this.img, param);
            
            % Update tracking points
            R = ECCWarp(:,1:2);
            t = ECCWarp(:,3);
            
%             for i=1:numel(this.segs)
%                 this.segs(i).transform('euclidean2D',[R,t]);
%             end
            this.segs.transform('euclidean2D',[R,t]);
        end
        function SVD_update( this )
            set( gcf, 'Name', 'SVD optimization' );
%             if this.WITH_PLOT
%                 clf
%                 imshow( this.img ); hold on;
%                 freezeColors;
%             end

            % Set mask of valid points (1 for valid)
            mask_segs = this.maskSegs;
            
            % Optimize for each segment
            Cpts  = cell(1,length(this.segs));
            Cm    = cell(1,length(this.segs));
            Cdist = cell(1,length(this.segs));
            newSegs = CSegment2D.empty(0,this.segs.Nsegs);
%             for i=1:numel(this.segs)
            for i=1:this.segs.Nsegs
                if ~mask_segs(i)
                    % If mask==0 for segment, skip update
                    continue
                end
                
                % Get the pixels in segment window
                % Window has width 2*h and length seg.d+2*g
                pts = this.segmentWindow( this.segs(i) );
                if isempty( pts )
                    warning('No points found in segment window');
                    keyboard;
                end
                
                % Optimize over obtained window
                [seg,pts] = this.optimSeg( this.segs(i), pts );
                
                % Get distance to center of points on new line
                m_ = [seg.v; -seg.v'*seg.c];
                dist_ = m_' * makehomogeneous( pts );
                if this.WITH_PLOT && this.WITH_DEB
                    % Plot old segment
%                     this.segs(i).plot('y',num2str(i));
                    % Plot new segment
                    seg.plot('c');
                end
                newSegs(i) = seg;
                
                % Store cell array of data for further refinement
                Cpts{i} = pts;
                Cm{i} = m_;
                Cdist{i} = dist_;
            end
            
            % Filter points in window using connected segments
            for i=1:this.segs.Nsegs
                if ~mask_segs(i)
                    % If mask==0 for segment, skip update
                    continue
                end
                [J1,J2] = this.segs.getNeighboursIdxs( i );
                % Remove masked segments from J1 and J2
                J1 = setdiff( J1, find(~mask_segs) );
                J2 = setdiff( J2, find(~mask_segs) );
                % Get intersection distance of current segment with
                % the rest of connected segments
                if ~isempty(J1)
                    dm = newSegs(i) * newSegs(J1) - newSegs(i);
                    maskm_ = Cdist{i} > max(dm);
                else
                    maskm_ = true( size(Cdist{i}) );
                end
                if ~isempty(J2)
                    dp = newSegs(i) * newSegs(J2) - newSegs(i);
                    maskp_ = Cdist{i} < min(dp);
                else
                    maskp_ = true( size(Cdist{i}) );
                end
                mask_ = maskm_ & maskp_;
                Cpts{i} = Cpts{i}(:,mask_);
            end
            
            % Global optimization applying constraints from connectivity
            % Optimize over cropped sets of points
            for i=1:this.segs.Nsegs
                if ~mask_segs(i)
                    % If mask==0 for segment, skip update
                    continue
                end
                [seg,pts] = this.optimSeg( newSegs(i), Cpts{i} );
                
                if this.WITH_PLOT && this.WITH_DEB
                    % Plot old segment
%                     this.segs(i).plot('y',num2str(i));
                    % Plot new segment
                    seg.plot('g');
                end
                newSegs(i) = seg;
            end
            
            % Compute points from multi-intersections
            for i=1:this.segs.Npts
                [J1,J2] = this.segs.getPointContainers( i );
                % Remove masked segments from J1 and J2
                J1 = setdiff( J1, find(~mask_segs) );
                J2 = setdiff( J2, find(~mask_segs) );
                if (numel(J1) + numel(J2)) < 2
                    % If there is no intersection, keep point
                    continue
                end
                
                segs1 = this.segs(J1);
                segs2 = this.segs(J2);
                L = [ segs1.l, segs2.l ];
                [~,~,V] = svd( L * L' );
                p_ = makeinhomogeneous( V(:,end) );
                this.segs.pts(:,i) = p_;
            end
            % Update object
            this.segs.updateSegs;
            
            % Remove previous graphic objects
            delete( this.hSegs ); this.hSegs = [];
            delete( this.hTags ); this.hTags = [];
            
            % Plot final points
%             plot(this.segs.pts(1,:), this.segs.pts(2,:), '*r');
            [this.hSegs,this.hTags] = this.segs.segs(mask_segs).plot('r',...
                this.segs.tags(mask_segs));
        end
        
        % Functions for optimization of a single segment
        function [seg, pts] = optimSeg( this, seg, pts )
            % Get indeces from array of points
            X = pts(1,:)'; Y = pts(2,:)';
            ind = sub2ind( this.imSize, Y(:), X(:) );
            mag = double( this.Gmag( ind ) );
            
            % Filter magnitude wrt median of all values
            mask_mag = this.filterMag(pts, median(mag));
            % Filter direction wrt 11.25º threshold from previous line
            mask_dir = this.filterDir(pts, seg.v, pi/16);
            % Fuse both masks
            mask = mask_mag & mask_dir;
%             mask = this.filterDir(pts, seg.v, pi/16);
            
            w = mag';
            l = svdAdjust( pts(:,mask), w(mask) );
                       
%             if this.WITH_DEB
%             imshow( this.img ); hold on;
%             freezeColors
%             this.plotSeg( seg );
%             end
            
            % Get projection of segment points on new line
            l1 = [ seg.v; -seg.v'*seg.p1 ];
            seg.p1 = makeinhomogeneous( cross( l, l1 ) );
            l2 = [ seg.v; -seg.v'*seg.p2 ];
            seg.p2 = makeinhomogeneous( cross( l, l2 ) );
            
%             if this.WITH_DEB
%             this.plotSeg( seg );
%             end

            % Check with DBSCAN
            % [C, ptsC, centres] = dbscan(dist, 5, 5);
            % dbscan too slow
            
            % Get distance of points on new line
            m = [seg.v; -seg.v'*seg.c];
            dist = m' * makehomogeneous( pts(:,mask) );
            
            prc1 = prctile( dist, 1 );
            prc2 = prctile( dist, 5 );
            prc3 = prctile( dist, 95 );
            prc4 = prctile( dist, 99 );
            c = seg.c;
            v = seg.v;
            seg.p1 = c + v * prc1; 
            seg.p2 = c + v * prc4;
            
%             if this.WITH_DEB
%             this.plotSeg( seg );
%             
%             figure, hist(dist), hold on
%             plot(prc1,0,'*r'), plot(prc2,0,'*r')
%             plot(prc3,0,'*r'), plot(prc4,0,'*r')
%             keyboard
%             end
            
            if this.WITH_DEB
%                 figure
%                 this.plotGdir( pts, seg );
%                 figure
                this.plotGdir( pts(:,mask), seg );
                this.plotSeg( seg );
                keyboard
                
%                 figure
%                 this.plotGmag( pts );
%                 figure
                this.plotGmag( pts(:,mask) );
                this.plotSeg( seg );
                keyboard
            end
            
            % Set output
            pts = pts(:,mask);
            
            function lin = svdAdjust( pts, w )
                % Computation (optimization) of line given by n weighted points in the segment
                % Input:
                %   pts - 2xN array of 2D points
                %   w - vector with weight of each point
                
                w = w(:);
                W = repmat( w', 3, 1 );
                
                pts = makehomogeneous( pts );
                
                pts = W .* pts;
                
                Q = pts * pts';
                
                Q_ = [ Q(1,1) - Q(1,3)^2 / Q(3,3)        , Q(1,2) - Q(1,3) * Q(2,3) / Q(3,3)
                       Q(1,2) - Q(1,3) * Q(2,3) / Q(3,3) , Q(2,2) - Q(2,3)^2 / Q(3,3)        ];
                [~,~,V] = svd( Q_ );
                lin = V(:,2);
                lin(3) = -(lin(1)*Q(1,3) + lin(2)*Q(2,3)) / Q(3,3);
            end
        end
        
        function [pts, corners] = segmentWindow( this, seg, h, g )
            % Get pixels inside orthogonal distance from segment
            
            % Assign inputs and variables
            p = seg.p1;
            q = seg.p2;
            n = seg.n; %#ok<*MCNPN>
            v = seg.v; %#ok<*MCNPN>
            % v is vector pointing from p to q
            
            imSize = this.imSize;
            if ~exist('h','var'), 
                h = this.h; %#ok<*PROP>
            end
            if ~exist('g','var')
                g = this.g; %#ok<*PROP>
            end
            
            % Find corner points for segment window
            corners = zeros(2,4);
            corners(:,1) = p + h * n - g * v;
            corners(:,2) = p - h * n - g * v;
            corners(:,3) = q + h * n + g * v;
            corners(:,4) = q - h * n + g * v;
            
            % Find container rectangle of current window (within image size)
            ymin = max( floor( min( corners(2,:) ) ), 1 );
            ymax = min(  ceil( max( corners(2,:) ) ), imSize(1) );
            xmin = max( floor( min( corners(1,:) ) ), 1 );
            xmax = min(  ceil( max( corners(1,:) ) ), imSize(2) );
            vx = xmin:xmax;
            vy = ymin:ymax;
            [X,Y] = meshgrid( vx, vy ); % All contained pixels indeces
            
            % Compute normalized homogeneous line
            l = lnormalize( cross( makehomogeneous(p), makehomogeneous(q) ) );
            
            % Get array of points
            pts = [ X(:) Y(:) ]';
            d   = ( l'  * makehomogeneous( pts ) );
            
            % Filter points according to orthogonal distance (constant distance)
            mask_d = abs(d) < h;
            pts = pts(:,mask_d);
        end
        function mask = filterMag( this, pts, thres )
            % mask = filterMag( this, pts )
            % pts is a 2xN array with pixel locations (indeces)
            % thres is the minimum gradient magnitude allowed
            % mask is a 1xN logical array 1 for pixels with ||G|| above
            % threshold
            
            % X is horizontal right (column, J index)
            % Y is vertical down (row, I index)
            X = pts(1,:)';
            Y = pts(2,:)';
            
            ind = sub2ind( this.imSize, Y, X );
            mask = this.Gmag(ind) >= thres;
            % Make column
            mask = mask(:);
        end
        function dtheta = angularDiff( this, pts, v )
            % dtheta = angularDiff( this, pts, v )
            % pts is a 2xN array with pixel locations (indexes)
            % v is a 2x1 S1 direction vector for image segment
            % This function computes the absolute angular difference
            % from 'unsigned' segment normal to grad directions
            % for each specified point in image
            
            % X is horizontal right (column, J index)
            % Y is vertical down (row, I index)
            X = pts(1,:)';
            Y = pts(2,:)';
            
            ind = sub2ind( this.imSize, Y, X );
            sindtheta = v' * this.Gdir(:,ind);
            dtheta = asin( abs(sindtheta) );
        end
        function mask = filterDir( this, pts, v, thres )
            % mask = filterMag( this, pts )
            % pts is a 2xN array with pixel locations (indexes)
            % v is a 2x1 S1 with segment direction
            % thres is the maximum angle difference permitted
            % mask is a 1xN logical array 1 for pixels with dtheta below
            % threshold
            
            % Obtain angular distance for gradient in points
            dtheta = this.angularDiff( pts, v );
            
            mask = dtheta <= thres;
            % Make column
            mask = mask(:);
        end
        
        % Debug functions
        function plotGdir( this, pts, seg )
            dtheta = this.angularDiff( pts, seg.v );
            this.plot( pts, dtheta );
            caxis([0 pi/4]);
            colorbar;
        end
        function plotGmag( this, pts )
            X = pts(1,:)'; Y = pts(2,:)';
            ind = sub2ind( this.imSize, Y(:), X(:) );
            mag = double( this.Gmag( ind ) );
            
            this.plot( pts, mag(:)' );
            caxis([0 max(mag)]);
            colorbar;
        end
        function plotSeg( this, seg )
            quiver( seg.p1(1), seg.p1(2),...
                    seg.v(1),  seg.v(2), seg.d );
        end
        
        function plot( this, pts, vals )
            figure(this.hFigure);
            imshow( this.img ); hold on;
            freezeColors
            
            scatter( pts(1,:), pts(2,:), 10, vals, 'filled' );
            colormap jet;
        end
        
        % Save and load methods
        function out = saveobj(this)
            
            % Make a struct with both interest properties (net and mask)
            out.segs = this.segs;
            out.maskSegs = this.maskSegs;
        end
    end
    methods (Static = true)
        function this = loadobj(data)
            this = CGtracker([],false);
            this.segs = data.segs;
            this.maskSegs = data.maskSegs;
        end
    end
  
end

