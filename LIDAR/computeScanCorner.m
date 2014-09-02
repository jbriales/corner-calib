function [l,A_l,A_lh,p,A_p,q,A_q, lin,seg] = computeScanCorner( scan, sigma, inliers, debug )
% [l,A_l,A_lh,p,A_p,q,A_q, lin,seg] = computeScanCorner( scan, inliers, debug )
% Input:
%   scan - scan structure with 2D points field (xy)
%   inliers - 1x3 cell array with inlier indexes for each segment
%   debug - show debugging
% 
% Output:
%   l,A_l - line direction and 1x1 angle covariance
%   A_lh  - homogeneous line covariance
%   p,A_p - 2D central point and 2x2 point covariance
%   q,A_q - 2D corner point
%   lin   - 1x3 cell array with homogeneous lines through segments
%   seg   - 1x3 cell array with 2 2x1 segment limit points

% Control variables
% sigma = 0.01; % SD of LIDAR sensor: Hokuyo UTM-30LX: 0.010 <10m and 0.030 <30m
% debug = 0;

% Input data
pts = scan.xy;

% Preallocation
[l,A_l,p,A_p,q,A_q,lin,seg] = deal( cell(1,3) );

%% Compute SVD adjustment and uncertainty propagation
% occ = cellfun(@(x)~isempty(x), {scantrack.lin});
occ = cellfun(@(x)~isempty(x), inliers);
for k=1:3
    if occ(k)
        seg_pts = pts(:,inliers{k});
        [l{k}, A_l{k}, p{k}, A_p{k}, lin{k}] = computeSegmentSVD( seg_pts, sigma, debug );
        
        % Compute segment end points
        seg{k}(:,1) = seg_pts(:,1) - lin{k}' * makehomogeneous( seg_pts(:,1) ) * lin{k}(1:2);
        seg{k}(:,2) = seg_pts(:,end) - lin{k}' * makehomogeneous( seg_pts(:,end) ) * lin{k}(1:2);
        
        % Compute segment biggest displacement from previous iteration
        % thres_disp(k) = max( abs( scantrack(k).lin' * makehomogeneous( out(k).seg ) ) ); % To apply yet
    else
        [l{k},A_l{k},p{k},A_p{k},lin{k},seg{k}] = deal([]);
    end
end

%% Compute calibration data from detected segments: segment directions (l), corner points (q)

% Compute homogeneous line uncertainty
for i = 1:3
    if ~isempty( p{i} ) && ~isempty( lin{i} ) 
        J_lh_l = [0 -1; 1 0; -p{i}(2) p{i}(1)];
        J_lh_p = [0 0 ; 0 0; lin{i}(2) lin{i}(1)];
        J_l = [0; 1];    
        A_l_ = J_l * A_l{i} * J_l';
        A_lh{i} = J_lh_l* A_l_ *J_lh_l' + J_lh_p* A_p{i} *J_lh_p';  
    else
        A_lh{i} = [];
    end
end

% Compute intersection points of three segments
% ----------------------------------------------
% Point intersection of lines in plane IJ and IK must lie in axis I (common
% axis)
occ = cellfun(@isempty, lin);
occ = ~occ;
for i=1:3
    idx = setdiff(1:3,i);
    if all( occ(idx) )
        q{i} = cross( lin{idx(1)}, lin{idx(2)} );
        q{i} = makeinhomogeneous( q{i} );
        % Compute corner points (q) uncertainty
        J_ = [-skew(lin{idx(2)}) skew(lin{idx(1)})];
        A_ = [A_lh{idx(2)} zeros(3,3);
              zeros(3,3) A_lh{idx(1)}];
        A_q{i} = J_ * A_ * J_';   
    else
        q{i} = [];
        A_q{i} = [];
    end
end

% debug = 0;
if debug
%     hF = figure('units','normalized','outerposition',[0 0 1 1]); hold on, axis equal;
    hF = figure( ); hold on, axis equal;
    plot(pts(1,:),pts(2,:),'.k')
    plotLIDARframe
    rgb = 'rgb';
    for k=1:3
        if ~isempty( lin{k} )
            plotHomLineWin( lin{k}, rgb(k) )
            plot(seg{k}(1,:),seg{k}(2,:),strcat('o',rgb(k)));
        end
    end

    for k=1:3
        if ~isempty( q{k} )
            plot(q{k}(1),q{k}(2), ['*',rgb(k)])
        end
    end
    
    axis equal
    keyboard
    close(hF)
end