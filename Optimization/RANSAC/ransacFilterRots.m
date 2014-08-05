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


function [R, inliers] = ransacFilterRots(linR, thres, feedback)
    
    if nargin == 2
	feedback = 0;
    end
    
    [rows, npts] = size(linR);
    
    if rows ~=9
        error('data is not 1x9 linearized rotation matrix');
    end
    
    if npts < 3
        warning('Only %d points',npts)
    end
    
    s = 1;  % Minimum No of points needed to define rotation.
        
    fittingfn = @(x)x;
    distfn    = @rotrotdist;
    degenfn   = @isdegenerate;

    [R, inliers] = ransac(linR, fittingfn, distfn, degenfn, s, thres, feedback);
    
    % Perform least squares fit to the inlying points
%     B = fitplane(linR(:,inliers));
    
%------------------------------------------------------------------------
% Function to define a rotation given its vectorized 9 components

% function R = id(X);
%     R = reshape(X,3,3);
    
%------------------------------------------------------------------------
% Function to calculate distances between a plane and a an array of points.
% The plane is defined by a 3x3 matrix, P.  The three columns of P defining
% three points that are within the plane.

function [inliers, R] = rotrotdist(R, X, t)
    
    npts = size(X,2);
    RR = repmat(R,1,npts);
    
    diff = RR-X;
    chordal_distance = sqrt( sum( diff.^2, 1 ) );
    d = 2 * asind( chordal_distance / (2*sqrt(2)) );
    
    inliers = find(abs(d) < t);
    
    
%------------------------------------------------------------------------
% Function to determine whether a set of 3 points are in a degenerate
% configuration for fitting a plane as required by RANSAC.  In this case
% they are degenerate if they are colinear.

function r = isdegenerate(X)
    r = false;