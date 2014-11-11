function [l, p, v, n] = adjustLine( pts )
% Function that solves direction of line

N = size(pts,2);

% Compute central point of segment
p     = mean(pts,2);

% Compute covariance matrix
pts_p = pts - repmat( p, 1, N );
M = pts_p * pts_p';

[~, ~, V] = svd(M);

% Get output values
v = V(:,1); % Line direction
n = V(:,2); % Normal direction
l = [ n; -n'*p ]; % Homogeneous line
end