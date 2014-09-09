function [R, t] = generate_random_poses( )
% Input options:
% N     - number of frames
% min_d - min distance to origin at which Camera is placed
% max_d - max distance to origin at which Camera is placed
% min_ang - min angular distance wrt canonical planes at which Camera is
% placed
% ang_z_max - max rotation wrt Z axis of camera pointing to origin
% ang_x_max - max rotation wrt X axis of camera pointint to origin

N = 1000;
min_d = 2;
max_d = 5;
min_ang = 10;
max_ang = 90 - min_ang;
ang_z_max = 2*pi;
ang_x_max = deg2rad( 30 );

% Function to normalize array of column vectors
normalize = @(X) X ./ repmat( sqrt( sum( X.^2, 1 ) ),3,1 );

%% World to Cam translation
% Directions from W to LIDAR
% Generate N valid values
cont = 0;
v_w_c = zeros(3,N);
t = zeros(3,N);
while cont < N
    v_w_c_ = rand( 3, 1 );
    v_w_c_ = normalize( v_w_c_ );
    % Distance from W to LIDAR
    d_w_c_  = min_d + (max_d-min_d) * rand( 1 );
    ang_ = acosd( v_w_c_ );
    if ~any( ang_ > max_ang )
        cont = cont + 1;
        v_w_c(:,cont) = v_w_c_;
        t(:,cont) = v_w_c_ * d_w_c_;
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
ang_z = ang_z_max * rand(1,N); % ang_z_max typically 2*pi
ang_x = ang_x_max * rand(1,N); % ang_x_max typically FOV
R = reshape( [R_x;R_y;R_z], 3,3,N );
for i=1:N
    % Check premultiplication and good order
    R(:,:,i) = R(:,:,i) * RotationX(ang_x(i)) * RotationZ(ang_z(i));
end

% Transform output to cell arrays:
t = mat2cell( t, 3, ones(1,N) );
R = mat2cell( R, 3, 3, ones(1,N) );
R = permute( R, [1 3 2] );

end