% RANSACFITPLANE - fits plane to 3D array of points using RANSAC
%
% Usage  [B, P, inliers] = ransacfitplane(XYZ, t, feedback)
%
% This function uses the RANSAC algorithm to robustly fit a plane
% to a set of 3D data points.
%
% Arguments:
%          XYZ - 3xNpts array of xyz coordinates to fit plane to.
%          t   - The distance threshold between data point and the plane
%                used to decide whether a point is an inlier or not.
%          feedback - Optional flag 0 or 1 to turn on RANSAC feedback
%                     information.
%
% Returns:
%           B - 4x1 array of plane coefficients in the form
%               b(1)*X + b(2)*Y +b(3)*Z + b(4) = 0
%               The magnitude of B is 1.
%               This plane is obtained by a least squares fit to all the
%               points that were considered to be inliers, hence this
%               plane will be slightly different to that defined by P below.
%           P - The three points in the data set that were found to
%               define a plane having the most number of inliers.
%               The three columns of P defining the three points.
%           inliers - The indices of the points that were considered
%                     inliers to the fitted plane.
%
% See also:  RANSAC, FITPLANE

% Copyright (c) 2003-2008 Peter Kovesi
% School of Computer Science & Software Engineering
% The University of Western Australia
% http://www.csse.uwa.edu.au/
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

% June 2003 - Original version.
% Feb  2004 - Modified to use separate ransac function
% Aug  2005 - planeptdist modified to fit new ransac specification
% Dec  2008 - Much faster distance calculation in planeptdist (thanks to
%             Alastair Harrison) 


function [R, inliers] = ransacFitTransNormals(corresps, thres, feedback)
    
    if nargin == 2
        feedback = 0;
    end
    
    [rows, npts] = size(corresps);
    
    if rows ~=5
        error('data is not 1x5 stack correspondence');
    end
    
    if npts < 3
        warning('Only %d points, 3 or more needed',npts)
    end
    
    s = 3;  % Minimum No of points needed to define rotation.
        
    fittingfn = @rotadjust;
    distfn    = @correspdist;
    degenfn   = @isdegenerate;

    [R, inliers] = ransac(corresps, fittingfn, distfn, degenfn, s, thres, feedback);
    
    % Perform least squares fit to the inlying points
%     B = fitplane(linR(:,inliers));
    
%------------------------------------------------------------------------
% Function to define a rotation given its vectorized 9 components

function R = rotadjust(X)
    N = X(1:3,:);
    L = X(4:5,:);
    R0 = [ 0 -1  0
           0  0 -1
           1  0  0 ];
    
    Lev_Fun = @(R) Fun( R, N, L );
    [ R, err, errNorm, W ] = LM_Man_optim(Lev_Fun,R0,'space','SO(3)','debug',0, 'maxIters', 200);
    
function [residual,J] = Fun( R, N, L )
residual = dot( N, R(1:3,1:2) * L, 1 )';
J = cross( R(1:3,1:2) * L, N, 1 )';

 
%------------------------------------------------------------------------
% Function to calculate distances between a plane and a an array of points.
% The plane is defined by a 3x3 matrix, P.  The three columns of P defining
% three points that are within the plane.

function [inliers, R] = correspdist(R, X, t)
    
    all_n = X(1:3,:);
    all_l = X(4:5,:);
    
    d = dot( all_n, R(:,1:2)*all_l );
    
    inliers = find(abs(d) < t);
    
    
%------------------------------------------------------------------------
% Function to determine whether a set of 3 points are in a degenerate
% configuration for fitting a plane as required by RANSAC.  In this case
% they are degenerate if they are colinear.

function r = isdegenerate(X)
    r = abs( det(X(1:3,:)) ) < 0.2;