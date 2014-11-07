function [x, A_x, imgtrack, checkImage ] = optimize( obj )

% Corner tracking algorithm in the first step where we initialise the
% gradient threshold vector

% Inputs:
% x0        - Initial estimation [px py theta_1 theta_2 theta_3]
% Q0        - Initial end points to find the corner
% img_gray  - Gray image in the current iteration

% Outputs:
% x         - Estimation of the Y-junction [px py theta_1 theta_2 theta_3]
% Q         - Current end points to find the corner

if ~exist('debug','var')
    debug = 0;
end

% Temporal interface to old function corner_calib
% Extract tracking struct data
x0 = [obj.meta.xi.c.X
      obj.meta.xi.v(1).x
      obj.meta.xi.v(2).x
      obj.meta.xi.v(3).x];
% Q0 = imgtrack.Q;
if ~isempty(obj.meta.mag_th)
    mag_th = obj.meta.mag_th;
else
    mag_th = [0 0 0];
end
if ~isempty(obj.meta.ang_th)
    ang_th = obj.meta.ang_th;
else
    ang_th = [0 0 0];
end

Q0 = cell(1,3);
for k=1:3
    if obj.expandLines(k)
        Q0{k} = [];
    else
        Q0{k} = obj.meta.q{k};
    end
end
% if obj.expandLines % Expand inlier points search to image borders
%     Q0 = [];
% else
%     Q0 = cell2mat( obj.meta.q ); % Use previous end points as limits for inlier points search
% end

% Parameters
h = 2;                    % width of the rectangles (all equals)
R_k = [h h h]';

dir_th = deg2rad(30.0);     % gradient direction threshold
lambda_ini = obj.lambda_ini;            % lambda min (distance to center in line direction)
lambda_end = 2000;          % lambda max (distance to center in line direction)
% TODO: make +inf or remove condition

p0 = x0(1:2);
p0 = p0(:); % Assure column vector

% [img_grad_mag, img_grad_dir] = imgradient(img_gray);
img_grad_mag = obj.Gmag;
img_grad_dir = obj.Gang;
% H = fspecial('gaussian',[5, 5], 0.5);
% img_gray_filtered = imfilter(img_gray, H);
% [~, img_grad_dir_filtered] = imgradient(img_gray_filtered);

% Preallocate lines in corner
tic
check_repeat = true; % TODO: Change loop exit condition
pts0 = {[],[],[]};
count = 1;

xi0 = Cxi( x0 ); obj.xi = Cxi( x0 );
while check_repeat && count < 10
    [x, pts, w, mask, ~, debug_out] = ...
        coarseTrack(obj, x0, lambda_ini, lambda_end, R_k,...
        img_grad_mag, img_grad_dir, mag_th, dir_th,...
        Q0, debug);
    if debug
        debugTracking(debug_out, img_gray)
    end
    % Store points in each iteration in object
    for i=1:3
        obj.pts{i} = pts{i}(:,mask{i}); % Store inlier points data in object for debug
        obj.w{i}   = w{i}(:,mask{i});   % Store weights for points in object for debug
    end
    % For debug:
    % hInlPts = obj.plotInlierPts
    % hInlW   = obj.plotWeightPts
    % hXiLines = obj.plotXi
    % Create xi object and compare distances
    obj.xi = Cxi( x );
    [d_c, d_ang] = obj.xi.distance( xi0 );
    check_repeat = d_c > 0.1 || any( d_ang > deg2rad(0.1) );
%     check_repeat = norm( x(1:2)-x0(1:2) ) > 0.1 || norm( x(3:5)-x0(3:5) ) > deg2rad( 0.1 ); % TODO: Improve condition
    x0 = x;
    xi0 = obj.xi;
    pts0 = pts;
    count = count + 1;
end

% TODO: Apply RANSAC here to impose line width of 2-3 pixels (all
% codirectional points out of this threshold should be rejected as an
% outlier)

% Update magnitude threshold
for i=1:3
%     mag_th(i) = median( w{i} ); % w is a cell array with grad mag values
%     mag_th(i) = prctile(w{i},80);
%     mag_th(i) = prctile(w{i},50);
    mag_th(i) = prctile(w{i},40);
%     mag_th(i) = prctile(w{i},75);
end
% Median is applied to the whole set of geometric area covered by trapezoid in image

% Estimate the solution with a Levenberg - Marquardt optimisation
for i=1:3
    pts{i} = pts{i}(:,mask{i});
    w{i}   = w{i}(:,mask{i});
end
obj.pts = pts; % Store inlier points data in object for debug
obj.w   = w;   % Store weights for points in object for debug

Fun = @(x) Corner_Fun( x, pts, w );
% [ x, err ] = LM_Man_optim(Fun, x0, 'space', 'Rn', 'weighted', true, 'debug',0);
[ xi, err ] = LM_Man_optim(Fun, xi0, 'space', 'Manifold',...
                           'weighted', true, 'debug',0);
obj.xi = xi;
% obj.xi = Cxi( x );

% Final xi:
% For debug:
% hInlPts = obj.plotInlierPts
% hInlW   = obj.plotWeightPts
% hXiLines = obj.plotXi

% Covariance estimation (TODO: With generic LSQ cov propagation)
[~, J, W] = Corner_Fun( xi, pts, w );
% It's difficult to take a source for the error when optimizing
% the segment. Since the optimization here is mainly done
% based on gradient magnitude, the noise should be applied
% to this variable.
% We need to stablish thus a photometric covariance for this value,
% and decide here to try 
% Old way, non-sense in using current W as covariance matrix!?
% i_W = diag( 1./diag(W) );
% A_x1 = inv( J' * i_W * J );

% Current proposal: Set noise in vector error function
% In this case F is vector of signed distances
% It continues to make non-sense, but it's temporal
% The noise should come propagated from pixel uncertainty (position)
% and gradient uncertainty (value, propagable from intensity values?)
A_g = eye( size(W) ); % Set non-correlated identity covariance matrix in distances (1 pixel typical deviation)
J_xi_F = (J'*W*J) \ J'*W;
A_x = J_xi_F * A_g * J_xi_F';

% Set output tracking structure
imgtrack.x = x;
imgtrack.mag_th = mag_th;
imgtrack.ang_th = dir_th;

% Find points further to center and take as new q points
c = x(1:2);
Q = cell(1,3);
for k=1:3
    Npts = size(pts{k},2);
    distance = sqrt( sum( (pts{k}-repmat(c,1,Npts)).^2, 1 ) );
    [max_dist, idx] = max(distance);
%     Q{k} = pts{k}(:,idx);
    Q{k} = c + obj.xi.v(k).X * max_dist;
end
imgtrack.q = cell2mat(Q);

%% Check line correctness with gradient direction
% Check in complete rectangle around detected segment
cov_angle = zeros(1,3);
for i=1:3
    ang = x(2+i);
    p = x(1:2);
    v = [ cos(ang) sin(ang) ]';
    n = [ -v(2) v(1) ]';
    h = R_k(i);
    
    q = getBorderIntersection( p, v, obj.img_size );
    [pts_rectangle, ~] = findClosePoints( p, q, h, n, lambda_ini, obj.img_size);

    X = pts_rectangle(1,:)';
    Y = pts_rectangle(2,:)';
    ind = sub2ind( obj.img_size, Y(:), X(:) );
    
    % Gradient angle filtering (on mag-filtered points):
    dir = double( img_grad_dir( ind ) );
    dir = deg2rad(dir);
    grad_dir = [ -cos(dir) , +sin(dir) ]'; % IMP: why cos need to be negative?
    R_im_dir = [ [0 1 ; -1 0]*v v ]; % Convert direction vectors to SR aligned with n
    grad_dir_ = R_im_dir' * grad_dir;
    th_ang = atan( grad_dir_(2,:)./grad_dir_(1,:) );
    
    cov_angle(i) = cov( th_ang ); % Or use abs?
end

% Auxiliar function
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

checkImage = false;

% if any(cov_angle > 0.10)
% if any(cov_angle > pi/8)
% if any(cov_angle > 0.750)
if 0
    if any(cov_angle > 0.90)
        warning('Check image segment angle')
        checkImage = true;
        keyboard
        
        J = [];
        Q = [];
    end
end

% Check segments are not over the same line
ang_dist_thres = deg2rad(2);
if abs(x(3)-x(4)) < ang_dist_thres ||...
   abs(x(4)-x(5)) < ang_dist_thres ||...
   abs(x(5)-x(3)) < ang_dist_thres
   warning('Are two lines the same?');
   checkImage = true;
end
end

function [E, J, W, U] = Corner_Fun( xi, Cpts, Cw )
% Extract elements of product manifold
c = xi.c.X; % Value of c in xi
v{1} = xi.v(1).X;
v{2} = xi.v(2).X;
v{3} = xi.v(3).X;

R90 = [ 0 -1 ;
        1  0 ];
for k=1:3
    % Extract block auxiliary variables
    pts = Cpts{k};
    Npts = size(pts,2);
    w   = Cw{k};
    n   = R90 * v{k}; % Normal vector
    cc  = repmat(c,1,Npts);
    
    e = ( n' * (pts - cc) )'; % Converted to column vector
    Jc   = - repmat( n, 1, Npts )';
    Jang = - (pts - cc)' * v{k};
    
    % Store blocks in cells
    Cell_E{k} = e'; % Stored as row for cell2mat later operation
    Cell_Jc{k}   = Jc';
    Cell_Jang{k} = Jang';
    Cell_w{k} = w; % Weights: Gradient magnitude
end
% The data was stored in transposed order
% to ease concatenation with cell2mat
E = cell2mat( Cell_E )';

J = cell2mat(Cell_Jc);
J = [ J
    blkdiag( Cell_Jang{:} ) ];
J = J';

W = diag( cell2mat( Cell_w ) );
end

function [E, J, W, U] = Corner_Fun_OLD( x, Cpts, Cw )

p_0 = x(1:2);

% N_i = [ size(pts{1},2), size(pts{2},2), size(pts{3},2) ];
% N   = N_i(end);
% E = zeros(N,1);
% J = zeros(5,N);
% W = zeros(N,N);
for i=1:3
    pts = Cpts{i};
    w   = Cw{i};
    ang = x(2+i);
    v = [+cos(ang), +sin(ang)]';
    n = [-sin(ang), +cos(ang)]';
    l = [ n', -n'*p_0 ]';
    d = ( l' * makehomogeneous( pts ) )';

    Npts = length(d);
    Cell_E{i} = d';
    Cell_J12{i} = repmat( -n, 1, Npts ); % Derivative wrt central point
%     Cell_J3{i}  = -(v_k(i,1)*(pts(1,mask)-p_0(1))+v_k(i,2)*(pts(2,mask)-p_0(2))); % Derivative wrt angle
%     Cell_J3{i}  = - ( cos(ang)*(pts(1,mask)-p_0(1)) + sin(ang)*(pts(2,mask)-p_0(2))); % Derivative wrt angle
    Cell_J3{i}  = - v' * ( pts - repmat( p_0,1,Npts ) ); % Derivative wrt angle
    Cell_w{i} = w; % Weights: Gradient magnitude
    
    % Parameters for uncertainty propagation
    C_xy(1,i+2)   = -cos(ang)*Npts;
    C_xy(2,i+2)   = +sin(ang)*Npts;
    C_yy(i+2,i+2) = sum( cos(ang)*(pts(2,:)-p_0(2))+ sin(ang)*(pts(1,:)-p_0(1)) );
end
    E = cell2mat( Cell_E )';
    
    J = cell2mat(Cell_J12);
    J = [ J
          blkdiag( Cell_J3{:} ) ];
    J = J'; % Need to transpose: der of k wrt x parameters
      
    W = diag( cell2mat( Cell_w ) );
    
    % U is the Jacobian of the implicit function [ x = f(pts,w) ]
    C_yy(1:2,3:5) = -C_xy(1:2,3:5);
    C_yy(3:5,1)   = +C_yy(1,3:5)';
    C_yy(3:5,2)   = -C_yy(2,3:5)';    
    
    C_xy = 2*C_xy;
    C_yy = 2*C_yy;
    
    U = -inv(C_yy)*C_xy';
end

function [x, Cell_pts, Cell_w, Cell_mask, mag_th, debug_out] = ...
    coarseTrack(obj, x, mu_ini, mu_end, R_k,...
    img_grad_mag, img_grad_dir, mag_th, dir_th,...
    Q, debug)

% Vectors are taken as columns
% Array of vectors are tiled columns
if ~exist('debug','var')
    debug = 0;
end

px = x(1);
py = x(2);
p_0 = [px; py];

size_img = size( img_grad_mag );

v = [ cos( x(3:5) ), sin( x(3:5) ) ]';
ort = [ 0 -1
        1  0 ];
n = ort * v;

v_k = v;
n_k = n;

for i = 1:3
    % Set current iteration values
    h = R_k(i);     % distance threshold (in pixels)
    v = v_k(:,i);
    n = n_k(:,i);   % normal vector to segment
%     q = N_k(:,i);   % is chosen at the border of image below
    ang = x(2+i);   % segment angle

    % Get q as intersection of current line with image borders
%     if isempty(Q)
    if isempty(Q{i})
        q = getBorderIntersection( p_0, v, size_img );
    else
%         q = Q(:,i); % increment length?
        q = Q{i}; % increment length?
    end

    kt = 1.5; % Growth rate in area width
    cont_growth = 1;
    Npts = 0;
    while Npts < 10 && cont_growth < 20
        % Get points within trapezoid with distances h1 and h2 from p and q
%         [pts, corners] = findClosePoints( p_0, q, [h 10*h], n, mu_ini, size_img); %#ok<NASGU>
        [pts, corners] = findClosePoints( p_0, q, [h h], n, mu_ini, size_img); %#ok<NASGU>
        % Security check
        if isempty(pts)
            warning('pts is empty in image tracking');
            [pts, corners] = findClosePoints( p_0, q, [h 10*h], n, mu_ini, size_img); %#ok<NASGU>
        end
        
        if 0
            mask_mag = obj.filterMag(pts, mag_th(i));
            mask_dir = obj.filterDir(pts, v, pi/16);
            mask = mask_mag & mask_dir;
        else
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
        end
        Npts = sum(mask);
        
        % Update parameters
        h0 = h;
        h  = kt * h;
        cont_growth = cont_growth + 1;
    end
    h = h0; % Recover last values of h
    
    % Adjust line to chosen points by SVD
    if 0
        
    else
        w = mag';
        l = svdAdjust( pts(:,mask), w(mask) );
    end
    
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

% Find intersection point of 3 lines
S = lin{1} * lin{1}' + lin{2} * lin{2}' + lin{3} * lin{3}';
p = - S(:,1:2) \ S(:,3); % Result of optimization

x = zeros(5,1);
x(1:2) = p;
ort = [  0 -1
        +1  0 ];
for k=1:3
    v = ort * lin{k}(1:2);
    c = mean(Cell_pts{k},2);
    pc = c - p;
    v = v * sign( v' * pc ); % Correct direction of v
    x(k+2) = atan2( v(2), v(1) );
end

if debug
    debug_out = {debug_pts, masks, img_grad_mag, debug_l1, debug_l2, debug_dir};
else
    debug_out = [];
end
end

% Auxiliar function for debugging
function debugTracking( Cdebug, img )
% Assign input arguments
Cpts = Cdebug{1};
masks = Cdebug{2};
% img = Cdebug{3};
l1 = Cdebug{4};
l2 = Cdebug{5};
grad_dir = Cdebug{6};

figure
% imshow( colorScale(img) ), hold on
imshow( img ), hold on

mask_mag = [ masks{1,1}', masks{2,1}', masks{3,1}' ];
% mask_diff = [ masks{1,2}', masks{2,2}', masks{3,2}' ];
% mask_lambda = [ masks{1,3}', masks{2,3}', masks{3,3}' ];
mask = [ masks{1,4}', masks{2,4}', masks{3,4}' ];

pts = cell2mat( Cpts );

% Different masks: Set from figure tools
h_mag = plot(pts(1,mask_mag),pts(2,mask_mag),['.','m']);
% h_diff = plot(pts(1,mask_diff),pts(2,mask_diff),['.','y']);
% h_lambda = plot(pts(1,mask_lambda),pts(2,mask_lambda),['.','c']);
h_d = plot(pts(1,:),pts(2,:),['.','w']);

set(h_mag,'Visible','off')
% set(h_diff,'Visible','off')
% set(h_lambda,'Visible','off')
set(h_d,'Visible','off')

%         plot(pts(1,mask),pts(2,mask),'.g')
rgb = 'rgb';
for i=1:3
    plot(Cpts{i}(1,masks{i,4}),Cpts{i}(2,masks{i,4}),[rgb(i),'.'])
    plotHomLineWin(l1{i},[rgb(i),'--'])
    plotHomLineWin(l2{i},[rgb(i),'-'])
    
    dir  = grad_dir{i};
    hdir(i) = quiver( Cpts{i}(1,masks{i,1}), Cpts{i}(2,masks{i,1}), dir(1,:), dir(2,:), 0.05, 'r' );
    n  = l2{i}(1:2);
    nn = repmat( n, 1, length(Cpts{i}(1,masks{i,4})) );
    hn(i) = quiver( Cpts{i}(1,masks{i,4}), Cpts{i}(2,masks{i,4}), nn(1,:), nn(2,:), 0.05, 'g' );
    hMask(i) = plot( Cpts{i}(1,masks{i,4}), Cpts{i}(2,masks{i,4}), 'yo' );
end
hLeg = legend('Grad mag', 'Grad dir', 'Distance to corner', 'Geometric distance',...
              'X points', 'X1', 'X2', 'Y points', 'Y1', 'Y2', 'Z points', 'Z1', 'Z2');
set(hLeg,'Visible','off')
end

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
    
function [pts, corners] = findClosePoints( p, q, h, n, mu, size_img)
% Get pixels inside orthogonal distance from segment

if numel(h)==1
    hp = h;
    hq = h;
elseif numel(h)==2
    hp = h(1);
    hq = h(2);
else
    error('Wrong nr of values in h vector')
end

% Find corner points
corners = zeros(2,4);
corners(:,1) = p + hp * n;
corners(:,2) = p - hp * n;
corners(:,3) = q + hq * n;
corners(:,4) = q - hq * n;

% Find container rectangle (within image size)
ymin = max( floor( min( corners(2,:) ) ), 1 );
ymax = min(  ceil( max( corners(2,:) ) ), size_img(1) );
xmin = max( floor( min( corners(1,:) ) ), 1 );
xmax = min(  ceil( max( corners(1,:) ) ), size_img(2) );
vx = xmin:xmax;
vy = ymin:ymax;
[X,Y] = meshgrid( vx, vy );

% Compute homogeneous line
l = cross( makehomogeneous(p), makehomogeneous(q) );
beta = sqrt(l(1)^2+l(2)^2);
l = l / beta; % Normalize line

% Get array of points
pts = [ X(:) Y(:) ]';
d   = ( l'  * makehomogeneous( pts ) );

if numel(h)==1
    % Filter points according to orthogonal distance (constant distance)
    mask_d = abs(d) < h;
    pts = pts(:,mask_d);
else
    % Filter points according to orthogonal distance (variable distance)
    t_max = norm(q-p);
    v = (q-p) / t_max;
    Npts = size(pts,2);
    t = v' * ( pts - repmat(p,1,Npts) );
    d_max = ( t*hq + (t_max-t)*hp ) ./ t_max;
    mask_d = abs(d) < d_max;
    pts = pts(:,mask_d);
end

% Distance to point p in line direction: lambda
v = q-p; v = v/norm(v);
% Use positive direction from p_0 only
N = size(pts,2);
pts_p = pts - repmat(p,1,N);
lambda = v' * pts_p;
% Filter points according to distance to p
mask_lambda = lambda > mu; % Take points only in positive direction from p

% Filter points that are further than original distance pq
if 1
%     warning('r distance filter is active');
    mask_r = lambda < norm(q-p) ;
else
    mask_r = ones(size(mask_lambda));
end
pts = pts(:,mask_lambda & mask_r);
end

function q = getBorderIntersection( p, v, size_img )
% Get intersection with image border of line going from p in direction v

n = [ -v(2) v(1) ]';
l = [ n' -n'*p ]';
ax = [ 1 size_img(2) 1 size_img(1) ];
borders = [ 1 1 0 0
    0 0 1 1
    -ax   ]; % Homogeneous lines for image borders
q = makeinhomogeneous( skew(l) * borders );
pq = q - repmat(p,1,4);
% Find closest q of the 4 results in positive v direction
q  =  q(:, v'*pq > 0);
pq = pq(:, v'*pq > 0 );  % Remove back intersections
[~,Iq] = min( v'*pq ); % Take closest intersection
q = q(:,Iq);
clear Iq
end