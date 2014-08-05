function [N, c, A_co, L_P2] = imgParams2nonCalibratedData( x )
% [N, c, A_co, hL] = pts2calibratedData( P, K )
% Input:
%   x - 1x5 array with: c point (in pixels), 3 angles of direction
%   K - intrinsic calibration matrix 
% Output:
%   N  - 2x3 array with direction of lines
%   c  - 2D corner point
%   A_co - covariance matrix with uncertainty of output data
%   L_P2 - 3x3 array whose columns are (hom.) calibrated lines

%% Set input values
c = x(1:2);
a(1) = x(3);
a(2) = x(4);
a(3) = x(5);

%% Compute homogeneous lines
l = cell(1,3);
for k=1:3
    n = [-sin(a(k)) cos(a(k))]';
    d = - n' * c;
    l{k} = [ n ; d ];
end

%% Output data
% Calibrated homogeneous lines
L_P2 = cell2mat( l );

% Minimal data for rotation
N = L_P2(1:2,:);

% Set uncertainty:
% TODO: make numerical estimation
% A_co = diag([ 0.01 0.01 0.01, 0.01 0.01 ].^2 );
A_co = [];

end