classdef CGoptimizer < handle
    %CGoptimizer Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        img
        Gmag    % Image gradient magnitude
        Gang    % Image gradient angle
        Gdir    % Image gradient direction
        img_size
        
        meta
        xi      % Object to be optimized
        
        % Parameters
        lambda_ini
        h   % Segment window size
        
        % Control parameters and behavior
        expandLines
        
        % Debug data
        pts
        w
    end
    
    methods
        function obj = CGoptimizer( img, meta )
            if ndims( img ) == 3 % Color image
                img = rgb2gray( img );
            end
            obj.img = img;
            obj.meta = meta;
            
            % Compute image gradients once
            [obj.Gmag,obj.Gang] = imgradient( img );
            ang = deg2rad( obj.Gang(:) );
            % Gdir is stored with linear indexing to avoid 3rd dimension
            obj.Gdir = [ -cos(ang) , +sin(ang) ]'; % Check sign! Image coordinates?
            obj.img_size = size( img );
            
            % Select thresholds
            obj.lambda_ini = 10;
            
            % Choose behavior with boolean properties
            obj.expandLines = [true, true, true]; % Search until image borders by default
            
%             obj.xi = meta.xi;
        end
        
        function [obj_xi, meta] = compute( obj )
            [xi, A_xi, imgtrack, CHECK_IMAGE] = obj.optimize;
            
            % Bridge to meta
            c = xi(1:2);
            q = num2cell(imgtrack.q,1);
            meta = CMetaImg( c, q, imgtrack.mag_th, imgtrack.ang_th );
            
            if CHECK_IMAGE
                obj_xi = [];
                meta.optimized = false;
                return
            end
            % If everything went well
            meta.optimized = true;
            obj_xi = Cxi( xi(1:2), xi(3), xi(4), xi(5) );
            obj_xi.setMinimalCov( A_xi ); % MinimalCov is 5x5 (2+1+1+1)
            meta.xi = obj_xi; % To store with covariance
        end
        
        [x, A_x, imgtrack, checkImage] = optimize( obj )
        
        function Gcoarse( obj, seg )
            % [x, Cell_pts, Cell_w, Cell_mask, mag_th, debug_out] = coarseTrack(x, mu_ini, mu_end, R_k, img_grad_mag, img_grad_dir, mag_th, dir_th, Q, debug)
            size_img = size( img_grad_mag );
            
            ort = [ 0 -1
                1  0 ];
            n = ort * v;
            
            % Get q as intersection of current line with image borders
            if isempty(Q)
                q = getBorderIntersection( p_0, v, size_img );
            else
                
            end
            
            kt = 1.5; % Growth rate in area width
            count_growth = 1;
            Npts = 0;
            while Npts < 10 && count_growth < 20
                % Get points within trapezoid with distances h1 and h2 from p and q
                [pts, corners] = findClosePoints( p_0, q, [h 10*h], n, mu_ini, size_img); %#ok<NASGU>
                % Security check
                if isempty(pts)
                    warning('pts is empty in image tracking');
                    [pts, corners] = findClosePoints( p_0, q, [h 10*h], n, mu_ini, size_img); %#ok<NASGU>
                end
                % TODO: look for best proportional constant
                %         d = ( l' * makehomogeneous( pts ) )';
                X = pts(1,:)';
                Y = pts(2,:)';
                
                % Extract gradient data
                ind = sub2ind( size_img, Y(:), X(:) );
                % Gradient magnitude filtering (on trapezoidal window points):
                mag = double( img_grad_mag( ind ) );
                mask_mag = mag > mag_th(i);
                % Gradient angle filtering (on mag-filtered points):
                dir = double( img_grad_dir( ind( mask_mag ) ) );
                dir = deg2rad(dir);
                grad_dir = [ -cos(dir) , +sin(dir) ]'; % IMP: why cos need to be negative?
                R_im_dir = [ [0 1 ; -1 0]*v v ]; % Convert direction vectors to SR aligned with n
                grad_dir_ = R_im_dir' * grad_dir;
                th_ang = atan( grad_dir_(2,:)./grad_dir_(1,:) );
                med_ang = median( th_ang ); % Use median to obtain a characteristic value of gradient angle
                diff_ang = th_ang - med_ang;
                mask_ang = abs(diff_ang) < pi/16; % Filter directions farther than 11.25deg
                
                % Set final mask concatenating filters:
                idx_mag = find( mask_mag );
                idx_ang = idx_mag( mask_ang );
                mask = false(length(mask_mag),1);
                mask(idx_ang) = true;
                Npts = sum(mask);
                
                % Update parameters
                h0 = h;
                h  = kt * h;
                count_growth = count_growth + 1;
            end
            h = h0; % Recover last values of h
            
            % Adjust line to chosen points by SVD
            w = mag';
            l = svdAdjust( pts(:,mask), w(mask) );
            
            % Remove points farther than selected segment width (outliers)
            idx_prev = find( mask );
            d = l' * makehomogeneous( pts(:,mask) );
            
            mask_rect = abs(d) > h;
            mask( idx_prev( mask_rect ) ) = false; % Remove outliers from rectangle
            
            lin{i} = l;
            
            if debug % Store debug parameters
                masks{i,1} = mask_mag;
                %         masks{i,2} = mask_diff;
                %         masks{i,3} = mask_lambda;
                masks{i,4} = mask;
                debug_pts{i} = pts;
                debug_l2{i} = lin{i};
                debug_l1{i} = l;
                debug_dir{i} = grad_dir;
            end
            
            %     Cell_pts{i} = pts(:,mask);
            Cell_pts{i} = pts;
            Cell_w{i}   = w;
            Cell_dir{i} = dir;
            Cell_mask{i} = mask;
        end
        
        % Auxiliar functions for debug and plotting
        % Auxiliar function
        function h = plotInlierPts( obj )
            color = 'rgb';
            hold on
            for k=1:3
                h(k) = plot( obj.pts{k}(1,:), obj.pts{k}(2,:), ['.',color(k)] );
            end
        end
        
        function h = plotWeightPts( obj )
            hold on
            for k=1:3
                % X,Y vectors
                % C a Nx3 matrix with each row the color of corresponding
                % point XY
                w = obj.w{k};
                wmax = max( w );
                N = size( obj.pts{k}, 2 );
                C = zeros(N,3);
                C(:,k) = w / wmax;
                h(k) = scatter( obj.pts{k}(1,:), obj.pts{k}(2,:), [], C); %#ok<AGROW>
            end
        end
        function h = plotXi( obj )
            hold on
            col = 'rgb';
            h = zeros(1,3);
            for k=1:3
                h(k) = plotHomLineWin( obj.xi.l(k), col(k) );
            end
        end
        
        function plotGradPts( pts, n, dir_ang )
            % Function for plotting gradient vector field and line normal in
            % given points
            vdir = [ -cos(dir_ang) , +sin(dir_ang) ]';
            hp = debugPlotPts( pts, 'r.' ); %#ok<NASGU>
            hv = quiver( pts(1,:), pts(2,:), vdir(1,:), vdir(2,:), 0.05, 'w' ); %#ok<NASGU>
            N = length(pts);
            nn = repmat( n, 1, N );
            hn = quiver( pts(1,:), pts(2,:), nn(1,:), nn(2,:), 0.05, 'r' ); %#ok<NASGU>
        end
        
        %% Useful modular tools
        function pts  = findWindow( obj, c, q )
            % pts  = findWindow( obj, c, q )
            % c is the vertex point in image
            % q is the segment vector end point in image
            % pts is a 2xN array of points coordinates inside image
            % Some object-specific properties are used when computed this
            % window:
            % h - the window width
            % lambda - the min distance of points to vertex
            
            h = obj.h; %#ok<*PROP>
            
            % Recover vector data from extreme points c and q
            v = snormalize( q-c );
            n = [-v(2), +v(1)]';
            
            corners = zeros(2,4);
            corners(:,1) = c + h * n;
            corners(:,2) = c - h * n;
            corners(:,3) = q + h * n;
            corners(:,4) = q - h * n;
            
            % Find container rectangle (within image size)
            ymin = max( floor( min( corners(2,:) ) ), 1 );
            ymax = min(  ceil( max( corners(2,:) ) ), size_img(1) );
            xmin = max( floor( min( corners(1,:) ) ), 1 );
            xmax = min(  ceil( max( corners(1,:) ) ), size_img(2) );
            vx = xmin:xmax;
            vy = ymin:ymax;
            [X,Y] = meshgrid( vx, vy );
            pts = [ X(:) Y(:) ]';
            
            % Filter wrt points normal distance to segment line
            l = cross( makehomogeneous(c), makehomogeneous(q) );
            d = l' * makehomogeneous( pts );
            mask_d = abs(d) < h;
            pts = pts(:,mask_d);
        end
        function mask = filterMag( obj, pts, thres )
            % mask = filterMag( obj, pts )
            % pts is a 2xN array with pixel locations (indexes)
            % thres is the minimum gradient magnitude allowed
            % mask is a 1xN logical array 1 for pixels with ||G|| above
            % threshold
            
            % X is horizontal right (column, J index)
            % Y is vertical down (row, I index)
            X = pts(1,:)';
            Y = pts(2,:)';
            
            ind = sub2ind( obj.img_size, Y, X );
            mask = obj.Gmag(ind) > thres;
            % Make column
            mask = mask(:);
        end
        function dtheta = angularDiff( obj, pts, v )
            % dtheta = angularDiff( obj, pts, v )
            % pts is a 2xN array with pixel locations (indexes)
            % v is a 2x1 S1 direction vector for image segment
            % This function computes the absolute angular difference
            % from 'unsigned' segment normal to grad directions
            % for each specified point in image
            
            % X is horizontal right (column, J index)
            % Y is vertical down (row, I index)
            X = pts(1,:)';
            Y = pts(2,:)';
            
            ind = sub2ind( obj.img_size, Y, X );
            sindtheta = v' * obj.Gdir(:,ind);
            dtheta = asin( abs(sindtheta) );
        end
        function mask = filterDir( obj, pts, v, thres )
            % mask = filterMag( obj, pts )
            % pts is a 2xN array with pixel locations (indexes)
            % v is a 2x1 S1 with segment direction
            % thres is the maximum angle difference permitted
            % mask is a 1xN logical array 1 for pixels with dtheta below
            % threshold
            
            % Obtain angular distance for gradient in points
            dtheta = obj.angularDiff( pts, v );
            
            mask = dtheta < thres;
            % Make column
            mask = mask(:);
        end
    end
    
end

