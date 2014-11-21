function [R, t] = P3P_LRF( q, V )
% [R, t] = P3P_LRF( P )
% Input:
%   q - 2x3 array whose columns are 2D points in LIDAR plane contained
%       in trihedron lines 1,2,3 respectively
%   V - 3x3 array whose columns are 3D directions in Camera frame corresponding
%       to trihedron lines 1,2,3 respectively
% Output:
%   R - rotation matrix of LRF SR wrt Camera SR
%   t - translation vector from trihedron vertex to LRF origin

% Compute squared distance between LRF points pairs
% Order of distances vector L2 is on planes YZ, ZX and XY
L2 = zeros(3,1);
for k=1:3
    ij = setdiff(1:3,k); i = ij(1); j = ij(2);
    L2(k) = norm( q(:,i) - q(:,j) )^2;
end

% Solve multipliers for trihedron vectors
% To get 3D position of intersection points (up to sign)
M = [0 1 1
     1 0 1
     1 1 0];

if norm(V'*V-eye(3),'fro') < 1e-10 % If orthogonal trihedron (check threshold)
    d2 = M \ L2;
    d  = sqrt(d2); % Absolute value of solution (undetermined sign)
else
    % Applying cosines rule
    C = V'*V;
    c = [ C(2,3), C(3,1), C(1,2) ];
    C = diag(c);
    
    % A P3P problem arises:
    % xi^2 + xj^2 -2xixj cos(theta_ij) = L_ij^2
    % M * [x2 y2 z2]' -2*C*[yz,zx,xy]' == L2
    % TODO: Solve P3P problem
end

% Get 3D points (vectors, with origin in trihedron vertex)
% in camera frame (up to sign):
P  = V * diag(d);
% Get homogeneous 2D points in LRF plane
hp = makehomogeneous( q );

% Solve existing transformation M = [Rx Ry t]
% Eight possible solutions (sign permutations) exist
count = 0;
R = cell(1,8);
t = cell(1,8);
for s1 = [-1,+1]
    for s2 = [-1,+1]
        for s3 = [-1,+1]
            count = count + 1;
            S = diag([s1 s2 s3]);
            M = P * S / hp;
            Rx = M(:,1);
            Ry = M(:,2);
            Rz = cross(Rx,Ry);
            R{count} = [Rx Ry Rz];
            t{count} = M(:,3);
        end
    end
end

end