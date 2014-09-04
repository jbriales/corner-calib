function [img] = computeCameraCorner( img, debug )
% [img] = computeScanCorner( img, debug )
% Input:
%   img     - img structure
%   debug   - show debugging
% 
% Output:
%   img     - img structure with the additional fields
%       x   - 5xN array with: c point (in pixels), 3 angles of direction
%       A_x - 5x5xN covariance matrices with the uncertainty of x

% Control variables
sigma = 0.5;    % SD of Camera sensor

% Define the x and A_x structures
N   = size(img.pts,3);
x   = zeros(5,N);
A_x = zeros(5,5,N);

% Assignment of auxiliar variables
c   = reshape(img.pts(:,1,:),2,N);
p   = [reshape(img.pts(:,2,:),2,N) - c; 
       reshape(img.pts(:,3,:),2,N) - c; 
       reshape(img.pts(:,4,:),2,N) - c];

% Estimation of x
x(1:2,:) = c;
for i = 1:3
    x(i+2,:) = atan2(p(2*i,:),p(2*i-1,:));
end

% Estimation of A_x
A_x(1,1,:) = sigma;
A_x(2,2,:) = sigma;

A_pc = sigma * eye(4);
J    = zeros(3,4,N);
for i = 1:3
    % Set the num and den variables
    px = p(2*i-1,:);
    py = p(2*i,:);
    % Set the derivatives of atan2
    alpha = px.^2 + py.^2;
    d_atan_x = -py./alpha;
    d_atan_y =  px./alpha;
    % Calculate the Jacobian matrix
    J(i,1,:) =  d_atan_x;
    J(i,2,:) =  d_atan_y;
    J(i,3,:) = -d_atan_x;
    J(i,4,:) = -d_atan_y;
end
for i = 1:N
    A_x(3:5,3:5,i) = J(:,:,i) * A_pc * J(:,:,i)';
end

% Assignment of x and A_x to img structure
img.x   =   x;
img.A_x = A_x;