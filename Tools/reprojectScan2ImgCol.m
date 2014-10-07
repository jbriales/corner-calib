function reprojectScan2ImgCol( im, xy, K, T_c_s, col )
% reprojectScan2Img( im, xy, K, T_c_s )
%   im - RGB image
%   xy - 2xN array of 2D points in LIDAR plane
%   K  - camera intrinsic calibration matrix
%   T_c_s - [R|t] euclidean transformation (LIDAR seen from Camera)

[uv, d] = reproject( im, xy, K, T_c_s );

% Show image with superimposed LIDAR points
for k=1:size(uv,2)
    plot(uv(1,k),uv(2,k), 'Color',col, 'Marker','.',...
        'LineWidth',1, 'MarkerSize',5);
%     pause
end
end

function [uv, d] = reproject( im, xy, K, T_c_s )

% Set parameters
R = T_c_s(1:3,1:3);
t = T_c_s(1:3,4);

N = size( xy, 2 );
[res_v, res_u, ~] = size(im);

% Transform LIDAR points to Camera reference frame
pts = R(:,1:2) * xy + repmat( t, 1, N );

% Transform 3D points to pixels
uv = K * pts;
uv = makeinhomogeneous( uv );

% Find indexes in list of laser points for which point is projected to
% image
mask = find( uv(1,:) >= 0 & uv(1,:) < res_u & ...
             uv(2,:) >= 0 & uv(2,:) < res_v );
uv = uv(:, mask);
d  = xy(1, mask);

end