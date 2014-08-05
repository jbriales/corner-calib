function [R_w_s, t_w_s, Ms] = co2LidarPose( P, s )
% function [R_w_s, t_w_s] = co2LidarPose( P, s )
% Input:
%   P - 2x3 array where columns are x, y and z 2D points in LIDAR plane
%   s - 3x1 array where values indicate sign of triplet points in World frame

% Compute squared distance between LIDAR triplet points
% Order of distances vector L2 is on planes YZ, ZX and XY
L2 = zeros(3,1);
vk = 1:3;
for k=1:3
    ij = setdiff(vk,k);
    i = ij(1); j = ij(2);
    L2(k) = norm( P(:,i) - P(:,j) )^2;
end

% Solve triplet coordinates in World frame
M = [0 1 1
     1 0 1
     1 1 0];
p2 = M \ L2;
p  = s .* sqrt(p2);
clear M

% Reparameterization of LIDAR plane and SR
% See Zisserman Part0.3 and notes on 07/07/14
% x = p(1);
% y = p(2);
% z = p(3);
XP3 = makehomogeneous( diag(p) );
xP2 = makehomogeneous( P );
% planeP3 = null( XP3 );
% M = null( planeP3' );
% G = (xP2' \ (XP3'*M))';
Ms = XP3 / xP2;
Ms = makeinhomogeneous( Ms );

t_w_s = Ms(:,3);
% R_w_s_x = Ms(:,1) - t_w_s;
R_w_s_x = Ms(:,1);
R_w_s_x = R_w_s_x / norm(R_w_s_x);
% R_w_s_y = Ms(:,2) - t_w_s;
R_w_s_y = Ms(:,2);
R_w_s_y = R_w_s_y / norm(R_w_s_y);
R_w_s = [ R_w_s_x R_w_s_y cross(R_w_s_x,R_w_s_y) ];

end