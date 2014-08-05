function [img, scan, gt] = generate_random_poses( varargin )
% [R_w_c, t_w_c, img_pts, scan_pts, all_outlier_idxs] = generate_random_poses( )
% Input (optionals):
%   N                   - number of frames  
%   sd_pixel            - standard deviation
%   sd_laser
%   outlier_k
%
% Output:
%   img                     - structure formed by:
%       img_pts             - 2x4xN array of 2D image points (Corner, X,Y and Z
%                             intersections with image border)
%       K                   - calibration matrix
%   scan                    - structure formed by:
%       scan_pts            - 2xscan_NxN points intersected with world planes
%       scan_inliers        - scan_NxN range values
%   gt                      - structure formed by:
%       R_w_c               - 3x3xN array of Cam orientations wrt World
%       t_w_c               - 3xN array of Cam positions wrt World
%       R_w_s               - 3x3xN array of Lidar orientations wrt World
%       t_w_s               - 3xN array of Lidar positions wrt World
%       R_c_s               - 3x3 matrix of Lidar orientation wrt Camera
%       t_c_s               - 3x1 array of Lidar position wrt Camera
%       all_outlier_idxs    - vector with indexes of Cam or Lidar outliers

close all  

%% Set input data
N = 3;
withPlot = true;

% World data recording contrains
min_d = 2;
max_d = 5;
min_ang = 10;
max_ang = 90 - min_ang;

% Camera configuration data
K = [248.321289 0 319.416809
     0 248.321289 249.839676
     0 0 1 ];
size_img = [480 640];

% Scan configuration data
scan_N = 1081;
scan_FOV = deg2rad( 270.2 );
scan_d_min = 0.1;
scan_d_max = 20;

% Calibration rig data
R_c_s = [ 0 0 1 ; -1 0 0 ; 0 -1 0 ]' * RotationX( deg2rad( 7 ) );
t_c_s = [0 -0.5 0]';

gt.R_c_s = R_c_s;
gt.t_c_s = t_c_s;

% Noise addition
sd_pixel = 0.5;     %sd_pixel = 0;
sd_laser = 0.02;    %sd_laser = 0;    
outlier_k = 5; % Proportional value for outlier distribution

% Optional input arguments
allfields = {'N',           N;
             'sd_pixel',    sd_pixel;
             'sd_laser',    sd_laser;
             'outlier_k',   outlier_k;
             'with_plot',   false;
            };
opt = setVarargin(allfields, varargin{:}); 

vars = fieldnames( opt );
for i=1:length(vars)
    field = vars{i};
    eval([ field,'= opt.',field,';' ]);
end

% FOV = deg2rad(30);
FOV = 0.25 * min( [2*atan(K(1,3)/K(1,1)), 2*atan(K(2,3)/K(2,2))] );

% Function to normalize array of column vectors
normalize = @(X) X ./ repmat( sqrt( sum( X.^2, 1 ) ),3,1 );

img_outlier_idxs  = 1:5:N;
scan_outlier_idxs = 1:7:N;
all_outlier_idxs  = union( img_outlier_idxs, scan_outlier_idxs );
gt.all_outlier_idxs = all_outlier_idxs;


%% World to Cam translation
% Directions from W to LIDAR
% Generate N valid values
cont = 0;
v_w_c = zeros(3,N);
t_w_c = zeros(3,N);
while cont < N
    v_w_c_ = rand( 3, 1 );
    v_w_c_ = normalize( v_w_c_ );
    % Distance from W to LIDAR
    d_w_c_  = min_d + (max_d-min_d) * rand( 1 );
    ang_ = acosd( v_w_c_ );
    if ~any( ang_ > max_ang )
        cont = cont + 1;
        v_w_c(:,cont) = v_w_c_;
        t_w_c(:,cont) = v_w_c_ * d_w_c_;
    end
end

%% World to Cam rotation
% Set firstly camera pointing to World origin with X axis in XY world plane
% direction
R_z = -v_w_c;
R_x = -skew([0 0 1]') * R_z;
R_x = normalize( R_x );
R_y = cross(R_z,R_x,1);
% Add random rotation on Cam X and Z axis
ang_z = 2*pi* rand(1,N);
ang_x = FOV * rand(1,N);
R_w_c = reshape( [R_x;R_y;R_z], 3,3,N );
T_Cam = zeros(4,4, N);
T_Cam(4,4,:) = 1;
for i=1:N
    % Check premultiplication and good order
    R_w_c(:,:,i)     = R_w_c(:,:,i) * RotationX(ang_x(i)) * RotationZ(ang_z(i));
    T_Cam(1:3,1:3,i) = R_w_c(:,:,i);
    T_Cam(1:3,4,i)   = t_w_c(:,i);
end
gt.R_w_c = R_w_c;
gt.t_w_c = t_w_c;

%% Image points
% Set image limits
ax = [ 1 size_img(2) 1 size_img(1) ];
borders = [ 1 1 0 0
            0 0 1 1
              -ax   ]; % Homogeneous lines for image borders
t_c_w = zeros(3,N);
R_c_w = zeros(3,3,N);
lines = zeros(3,3,N); % Homogeneous lines for axis
img_q = zeros(2,3,N); %
img_pts = zeros(2,4,N);
for i=1:N
    t_c_w(:,i)      = -R_w_c(:,:,i)' * t_w_c(:,i);
    R_c_w(:,:,i)    =  R_w_c(:,:,i)';
    % Points are represented in pixel space
    c = K * t_c_w(:,i);
    inf_pts = K * R_c_w(:,:,i);
    lines_ = skew(c) * inf_pts;
    lines(:,:,i) = lines_;
    
    % Intersect lines with image borders to get observed points
    v = makeinhomogeneous(inf_pts) - repmat(makeinhomogeneous(c),1,3);
    for k=1:3
        %IMP: check sign of v(:,k)
        img_q(:,k,i) = getBorderIntersection( makeinhomogeneous(c), -v(:,k), lines_(:,k), borders );
    end
    img_pts(:,:,i) = [ makeinhomogeneous(c), img_q(:,:,i) ];
end

K_inv = inv(K);
for i = 1:N
    for k = 1:4
        img_pts_repr(:,k,i) = K_inv * 5*[img_pts(:,k,i); 1];   
        pts_repr(:,k,i)     = R_w_c(:,:,i) * img_pts_repr(:,k,i) + t_w_c(:,i);
    end    
end



% Corrupt with gaussian noise
img_pts = img_pts + sd_pixel * randn(2,4,N);
for i=img_outlier_idxs
    img_pts(:,:,i) = img_pts(:,:,i) + outlier_k * sd_pixel * (-1+2*rand(2,4)); % Uniform distribution for outlier
end
img.pts = img_pts;
img.K   = K;
 
%% LIDAR poses and plane in World coordinates
t_w_s = zeros(3,N);
R_w_s = zeros(3,3,N);
plane = zeros(4,N);
Ms    = zeros(4,3,N);
for i=1:N
    t_w_s(:,i) = t_w_c(:,i) + R_w_c(:,:,i) * t_c_s;
    R_w_s(:,:,i) = R_w_c(:,:,i) * R_c_s;
    n_plane = R_w_s(:,3,i);
    d_plane = n_plane' * t_w_s(:,i);
    plane(:,i) = [ n_plane ; -d_plane ]; % Plane in P3
end
gt.R_w_s = R_w_s;
gt.t_w_s = t_w_s;

%% LIDAR rays and collision
% scan_N = 30; % Debug
% scan_res = scan_FOV / scan_N;
scan_theta = linspace( -scan_FOV/2, +scan_FOV/2, scan_N );
dir = [ cos( scan_theta )' sin( scan_theta )' ]';
scan_lrays = [ -sin( scan_theta )' cos( scan_theta )' zeros(scan_N,1) ]';
% Get lines of intersection with World planes in parametrized LIDAR plane
plane_lines = zeros(3,3,N); % Each column is the intersection of LIDAR plane with k-th World plane
world_planes = [ eye(3) ; 0 0 0 ];
line_normalize = @(X) X ./ repmat( sqrt( sum( X(1:2,:).^2, 1 ) ),3,1 );
scan_r   = zeros(scan_N, N);
scan_pts = zeros(2,scan_N, N);
scan_inliers = cell(3,N);
D = zeros( 3, scan_N );
for i=1:N
    Ms(:,:,i) = makehomogeneous( [ R_w_s(:,1:2,i) t_w_s(:,i) ] ); % Parametrized plane matrix in P2
    plane_lines(:,:,i) = line_normalize( Ms(:,:,i)' * world_planes ); 

%     figure, hold on, rgb = 'rgb'; axis equal, xlabel('X'), ylabel('Y')
% 	  Os = zeros(2,scan_N);
%     quiver(Os(1,:),Os(2,:),dir(1,:),dir(2,:))
%     plot(0,0,'>r')
    for k=1:3 % Better parallelization
        x = skew( plane_lines(:,k,i) ) * scan_lrays;
        x = makeinhomogeneous( x );
        d = dot( dir, x, 1 );
        D(k,:) = d;
%         debugPlotPts( x, ['.',rgb(k)] )
%         plotHomLineWin( plane_lines(:,k,i), rgb(k) )
    end
    
    D( D<0 ) = inf;
    [scan_r_, scan_label] = min( D );
    scan_r_( scan_r_ < scan_d_min | scan_r_ > scan_d_max ) = 0; % Crop invalid values out of distance
    scan_label( scan_r_ < scan_d_min | scan_r_ > scan_d_max ) = 0; % Crop invalid values out of distance
    for k=1:3
        scan_inliers{k,i} = find( scan_label==k );
    end
    scan_r_ = scan_r_ + sd_laser * randn(1,scan_N);
    if ~isempty( i == img_outlier_idxs )
        scan_r_ = scan_r_ + outlier_k * sd_laser * (-1+2*rand(1,scan_N));
    end
    scan_r(:,i) = scan_r_;
    scan_pts(:,:,i) = repmat(scan_r(:,i)',2,1).*dir;
%     debugPlotPts( scan_pts(:,:,i), '.k' )
end
scan.pts     = scan_pts;
scan.inliers = scan_inliers;

% % Corrupt with gaussian noise
% scan_pts = scan_pts + sd_laser * randn(2,scan_N,N);
% for i=img_outlier_idxs
%     scan_pts(:,:,i) = scan_pts(:,:,i) + outlier_k * sd_laser * (-1+2*rand(2,scan_N)); % Uniform distribution for outlier
% end

% %% Intersection of LIDAR plane with axis: q points
% % All axis orthogonal
% scan_q = zeros(3,3,N);
% I = eye(4);
% for i=1:N
%     for k=1:3
%         ij = setdiff(1:3,k);
%         M = [ plane(:,i) I(:,ij(1)) I(:,ij(2)) ];
%         hom_q = null(M');
%         scan_q(:,k,i) = makeinhomogeneous( hom_q );
%     end
% end

% % Check solution with data here
% n = reshape(R_c_w,3,[]);
% l = reshape( plane_lines(1:2,:,:), 2,[] );
% ort = [0 -1 ; 1 0];
% l = ort*l;
% res = dot(n, R_c_s(:,1:2)*l, 1)
% 
% gt.R_c_w = R_c_w;
% gt.l = l;

%% Representation of results
if withPlot
    % Plot image points
    figure('Name', 'Projection of image points'), hold on
    for i=1:N
        scatter(img_pts(1,:,i),img_pts(2,:,i))
    end
    axis( ax );
    
    % Plot camera poses
    figure
    scatter3( t_w_c(1,:), t_w_c(2,:), t_w_c(3,:), 'r' )
    hold on
    plotframe(eye(4),5,'W','m')
    for i=1:N
        %     T = eye(4);
        %     T(1:3,1:3) = [R_x(:,i) R_y(:,i) R_z(:,i)];
        %     T(1:3,1:3) = R(:,:,i);
        %     T(1:3,4) = t_w_c(:,i);
        plotframe( T_Cam(:,:,i), 0.2, '', 'rgb' );
    end
    axis equal
    rotate3d on
    
    % Plot LIDAR poses
    scatter3( t_w_s(1,:), t_w_s(2,:), t_w_s(3,:), 'b' )
    hold on
    plotframe(eye(4),5,'W','m')
    for i=1:N
        T_scan = eye(4);
        T_scan(1:3,1:3) = R_w_s(:,:,i);
        T_scan(1:3,4) = t_w_s(:,i);
        plotframe( T_scan, 0.2, num2str(i), 'myc' );
        
        rig_line = [t_w_c(:,i) t_w_s(:,i)]';
        line( rig_line(:,1), rig_line(:,2), rig_line(:,3), 'Color', 'k', 'LineWidth', 2 )
    end
    axis equal
    rotate3d on
    
   % Plot LIDAR points
    for i=1:N
        pts = repmat( t_w_s(:,i),1,scan_N ) + R_w_s(:,1:2,i) * scan_pts(:,:,i);
        scatter3( pts(1,:), pts(2,:), pts(3,:), '.k' )
        rep_t = repmat(t_w_s(:,i),1,scan_N);
        dir3  = R_w_s(:,1:2,i) * dir;
        hQ = quiver3( rep_t(1,:), rep_t(2,:), rep_t(3,:),...
            dir3(1,:),dir3(2,:),dir3(3,:), 0.2 );
        keyboard
        delete( hQ )
        clear pts
    end
    
    % Plot CAMERA points
    for i=1:N
        
        pts_aux1 = pts_repr(:,1,i);
        t = t_w_c(:,i);
        for k = 2:4
            pts_aux = [pts_aux1 pts_repr(:,k,i)];
            line(pts_aux(1,:) , pts_aux(2,:) , pts_aux(3,:) );
            clear pts_aux;
        end
        col = 'krgb';
        for k = 1:4
            pts_aux = [t pts_repr(:,k,i)];
            line(pts_aux(1,:) , pts_aux(2,:) , pts_aux(3,:),'Color',col(k) );
            clear pts_aux;
        end
%         pts = repmat( t_w_s(:,i),1,scan_N ) + R_w_s(:,1:2,i) * scan_pts(:,:,i);
%         scatter3( pts(1,:), pts(2,:), pts(3,:), '.k' )
%         rep_t = repmat(t_w_s(:,i),1,scan_N);
%         dir3  = R_w_s(:,1:2,i) * dir;
%         hQ = quiver3( rep_t(1,:), rep_t(2,:), rep_t(3,:),...
%             dir3(1,:),dir3(2,:),dir3(3,:), 0.2 );
%         keyboard
%         delete( hQ )
%         clear pts
    end
end

end

function q = getBorderIntersection( c, v, l, borders )
% Get intersection with image border of line going from p in direction v
q = makeinhomogeneous( skew(l) * borders );
cq = q - repmat(c,1,4);
% Find closest q of the 4 results in positive v direction
q  =  q(:, v'*cq > 0);
cq = cq(:, v'*cq > 0 );  % Remove back intersections
[~,Iq] = min( v'*cq ); % Take closest intersection
q = q(:,Iq);
clear Iq
end

function h = plotframe(T, len, label, color)

    if ~all(size(T) == [4,4])
        if all(size(T) == [3,3])
            z = eye(4);
            z(1:3,1:3) = T;
            T = z;
        elseif all(size(T) == [3,4])
            z = eye(4);
            z(1:3,1:4) = T;
            T = z;
        else
            error('plotframe: matrix is not 4x4')
        end
    end
    
    if ~exist('len','var')
        len = 1;
    end
    
    if ~exist('label','var')    
        label = '';
    end
    
    if ~exist('color','var')    
        color = 'k';
    else
        if length(color)==1
            color = repmat(color,1,3);
        end 
    end    
    
    % Assume scale specified by T(4,4) == 1
    
    origin = T(1:3, 4);             % 1st three elements of 4th column
    X = origin + len*T(1:3, 1);     % point 'len' units out along x axis
    Y = origin + len*T(1:3, 2);     % point 'len' units out along y axis
    Z = origin + len*T(1:3, 3);     % point 'len' units out along z axis
    
    h(1) = line([origin(1),X(1)], [origin(2), X(2)], [origin(3), X(3)], 'color', color(1));
    h(2) = line([origin(1),Y(1)], [origin(2), Y(2)], [origin(3), Y(3)], 'color', color(2));
    h(3) = line([origin(1),Z(1)], [origin(2), Z(2)], [origin(3), Z(3)], 'color', color(3));
    h(4) = line([X(1),Y(1)], [X(2), Y(2)], [X(3), Y(3)], 'color', color(3));
    
    h(5) = text(X(1), X(2), X(3), ['x' label], 'color', color(1));
    h(6) = text(Y(1), Y(2), Y(3), ['y' label], 'color', color(2));
    h(7) = text(Z(1), Z(2), Z(3), ['z' label], 'color', color(3));
end
