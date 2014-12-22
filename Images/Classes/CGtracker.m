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
        
        % Tracking parameters
        h = 3
    end
    
    methods
        function this = CGtracker( img )
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
            load( file, '-mat' );
            this.segs = temp_segs;
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
                end
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
            
            for i=1:numel(this.segs)
                this.segs(i).transform('euclidean2D',[R,t]);
            end
        end
        function SVD_update( this )
            set( gcf, 'Name', 'SVD optimization' );
            for i=1:length(this.segs)
                this.segs(i) = this.optimSeg( this.segs(i) );
                % Plot new segment
                this.segs(i).plot('g',num2str(i));
            end            
        end
        
        % Functions for optimization of a single segment
        function seg = optimSeg( this, seg )
            % Get the pixels in segment window
            pts = this.segmentWindow( seg );
            if isempty( pts )
                warning('No points found in segment window');
                keyboard;
            end
            
            % Get indeces from array of points
            X = pts(1,:)'; Y = pts(2,:)';
            ind = sub2ind( this.imSize, Y(:), X(:) );
            mag = double( this.Gmag( ind ) );
            
%             mask_mag = this.filterMag(pts, mag_th(i));
%             mask_dir = this.filterDir(pts, v, pi/16);
%             mask = mask_mag & mask_dir;
            mask = this.filterDir(pts, seg.v, pi/16);
            
            w = mag';
            l = svdAdjust( pts(:,mask), w(mask) );
            
            % Get projection of segment points on new line
            l1 = [ seg.v; -seg.v'*seg.p1 ];
            seg.p1 = makeinhomogeneous( cross( l, l1 ) );
            l2 = [ seg.v; -seg.v'*seg.p2 ];
            seg.p2 = makeinhomogeneous( cross( l, l2 ) );
            
%             figure
%             this.plotGdir( pts, seg );
%             
%             figure
%             this.plotGdir( pts(:,mask), seg );
%             keyboard
            
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
        
        function [pts, corners] = segmentWindow( this, seg )
            % Get pixels inside orthogonal distance from segment
            
            % Assign inputs and variables
            p = seg.p1;
            q = seg.p2;
            n = seg.n; %#ok<*MCNPN>
            
            imSize = this.imSize;
            h = this.h; %#ok<*PROP>
            
            % Find corner points for segment window
            corners = zeros(2,4);
            corners(:,1) = p + h * n;
            corners(:,2) = p - h * n;
            corners(:,3) = q + h * n;
            corners(:,4) = q - h * n;
            
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
            
            ind = sub2ind( this.img_size, Y, X );
            mask = this.Gmag(ind) > thres;
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
            
            mask = dtheta < thres;
            % Make column
            mask = mask(:);
        end
        
        % Debug functions
        function plotGdir( this, pts, seg )
            
            dtheta = this.angularDiff( pts, seg.v );
            this.plot( pts, dtheta );
        end
        function plot( this, pts, vals )
            imshow( this.img ); hold on;
            freezeColors
            
            scatter( pts(1,:), pts(2,:), 10, vals, 'filled' );
            colormap jet;
            caxis([0 pi/4]);
            
            % Set image axes
%             axis([1 this.imSize(1) 1 this.imSize(2)]);
%             axis ij
        end
    end
  
end

