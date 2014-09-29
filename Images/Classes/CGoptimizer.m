classdef CGoptimizer
    %CGoptimizer Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        img
        Gmag    % Image gradient magnitude
        Gdir    % Image gradient direction
        
        meta
        xi      % Object to be optimized
    end
    
    methods
        function obj = CGoptimizer( img, meta )
            if ndims( img ) == 3 % Color image
                img = rgb2gray( img );
            end
            obj.img = img;
            obj.meta = meta;
            
%             [obj.Gmag,obj.Gdir] = imgradient( img );
            
%             obj.xi = meta.xi;
        end
        
        function [obj_xi, meta] = compute( obj )
            % Temporal interface to old function corner_calib
            meta = obj.meta;
            img  = obj.img;
            
            imgtrack.x = [meta.xi.c.X
                meta.xi.v(1).x
                meta.xi.v(2).x
                meta.xi.v(3).x];
            if ~isempty(meta.mag_th)
                imgtrack.mag_th = meta.mag_th;
            else
                imgtrack.mag_th = [0 0 0];
            end
            if ~isempty(meta.ang_th)
                imgtrack.ang_th = meta.ang_th;
            else
                imgtrack.ang_th = [0 0 0];
            end
            imgtrack.q = cell2mat( meta.q );
            
            debug = 0;
            [xi, A_xi, imgtrack, CHECK_IMAGE] = ...
                corner_calib(imgtrack, img, debug);
            
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
    end
    
end

