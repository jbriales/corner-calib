function [lin, inliers_all] = UpdateInliers(inliers, all_pts, thres, prev_line, debug)
% [lin, inliers] = UpdateInliers(inliers, pts, thres)
% Take previous inliers and LIDAR scan and refine list of inliers with RANSAC
%   inliers_all - list of inlier indexes in complete scan
%   inliers_reduced - list of inlier indexes in reduced set of points

N = length(inliers);
growth = 0.5;
margin = ceil( 0.5*growth*N );
inliers_all = (inliers(1) - 2*margin) : (inliers(end) + 2*margin);

% Crop indexes below 1 and over max (of size of X)
NX = size(all_pts,2);
inliers_all(1:find(inliers_all<1,1,'last')) = [];
inliers_all(find(inliers_all>NX,1,'first'):end) = [];

% Filter inliers whose distance is larger than expected step
[inliers_prev_line] = line2ptDist( prev_line, all_pts(:,inliers_all), 100*thres ); %TODO: Use previous iteration change as threshold
inliers_all = inliers_all( inliers_prev_line );

% Get refined line and inliers by RANSAC
% Use only points of incremented set of initial inliers (inliers_all)
pts = all_pts( :, inliers_all );
[lin, inliers_red] = ransacfitline2D(pts, thres, 0);

% Remove biased points of head and tail of inliers:
inliers_red = removeHeadTail( pts, inliers_red, lin, thres, debug );

% Final selection of inliers in complete scan
inliers_all = inliers_all(inliers_red); % Take good inliers of incremented initial set

end

function [inliers] = line2ptDist(l, X, t)
    l = l / norm(l(1:2));    
    d = abs( l' * makehomogeneous(X) );
    inliers = find(d < t);
end

function inliers = removeHeadTail( pts, inliers, lin, thres, debug )

% Remove biased points of head and tail of inliers:
% Debug:
% debug = 1;
if debug
    figure, hold on, axis equal
    plot( pts(1,inliers), pts(2,inliers), '.-k' );
    plotHomLineWin( lin, 'b' )
end

d = lin' * makehomogeneous( pts(:,inliers) );
s = sign(d);
s_head = sign(d(1));
s_tail = sign(d(end));
head = find( s ~= s_head, 1, 'first' );
tail = find( s ~= s_tail, 1, 'last' );

% Control cases without losing
if abs(d(1)) < 0.5 * thres
    head = 1;
end
if abs(d(end)) < 0.5 * thres
    inliers = inliers(head:end);
else
    inliers = inliers(head:tail);
end

% Debug:
if debug
    plot( pts(1,inliers), pts(2,inliers), 'or' );
end
end