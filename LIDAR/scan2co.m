function [l,A_l,A_lh,p,A_p,q,A_q, scantrack, seg_pts, lost] = scan2co( inputScan, scantrack, debug )
% function [l,A_l, p,A_p, q, lin, seg, inliers, outPts, lost] = scan2co( inputScan, scantrack, debug )
% [l,A_l, p,A_p, q, lin, outPts] = scan2co( inputScan, scantrack, debug )
% Output:
%   l - direction vector (S1)
%   A_l - covariance of directions
%   A_lh - covariance of homogeneous line
%   p - central point in segment (R2)
%   A_p - covariance of central points
%   q - corner points in scan
%   A_q - covariance of corner points in scan
%   lin - homogeneous lines (P2)
%   seg - extreme points defining scan (two 2D points)
%   outPts - set of inliers (2D points) in each segment
% Note: All outputs are cell arrays so that some of them can remain empty
% (in the general corner observation case)

if ~exist('debug','var')
    debug = false;
end

% Control variables
threshold = 0.02;
min_inliers  = 30;
min_inliers2 = 5; % For secundary pieces of line
sigma = 0.01; % SD of LIDAR sensor: Hokuyo UTM-30LX: 0.010 <10m and 0.030 <30m
% debug = 0;

% Input data
pts = inputScan.xy;
x = pts(1,:);
y = pts(2,:);

seg_pts = cell(1,3);
lost = false(1,3);
% new_inliers = cell(1,3);
lin = cell(1,3);

%% Update segments and inliers
% occ = cellfun(@(x)~isempty(x), inputScan.line);
occ = cellfun(@(x)~isempty(x), {scantrack.lin});
for k=1:3
    if occ(k)
%         inliers = inputScan.inliers{k};
        % Update each line piece separately
        for p = 1:scantrack(k).n
            inliers = scantrack(k).inliers{p};
            prev_line = scantrack(k).lin;
            [lin{k}, new_inliers] = UpdateInliers(inliers, pts, threshold, prev_line, debug);
            
            % TODO: Check that piece of line is not lost
            if length(new_inliers) < min_inliers
                if p==1
                    warning('SCAN2CO: Main piece of segment is lost')
                    scantrack(k).n = scantrack(k).n - 1; % Decrease number of pieces
                    new_inliers = [];
                elseif length(new_inliers) < min_inliers2
                    warning('SCAN2CO: One piece of segment is lost')
                    scantrack(k).n = scantrack(k).n - 1; % Decrease number of pieces
                    new_inliers = [];
                end
            end
            scantrack(k).inliers{p} = new_inliers;
        end
        
        % Old code:
        %         lin = inputScan.line{k};
        %         scan(k) = RefineSegments(lin,x,y,threshold,min_inliers,sigma,debug);
        %         seg = inputScan.seg{k};
        %         [scan(k), lost(k)] = RefineSegmentsInliers(inliers,x,y,threshold,min_inliers,sigma,debug);
        %         scan(k) = RefineSegments(seg,x,y,threshold,min_inliers,sigma,debug);
    end
end

%% Check that two segments are not becoming the same
% If so, set lost the smallest segment
for k=1:3
    ij = setdiff(1:3,k);
    i = ij(1); j = ij(2);
    if occ(i) && occ(j)
        if rad2deg( acos( lin{i}(1:2)' * lin{j}(1:2) ) ) < 15 % Threshold: TODO: Check best value
            lab = 'XYZ';
            warning('Segments %c and %c are too similar',lab(i),lab(j))
%             Ni = length( inputScan.inliers{i} );
%             Nj = length( inputScan.inliers{j} );
            Ni = scantrack(i).Ninliers;
            Nj = scantrack(j).Ninliers;
            if Ni<Nj
                ind = i;
            else
                ind = j;
            end
            occ(ind) = false;
            lin{ind} = [];
%             new_inliers{ind} = []; % Non sense with multisegment
            lost(ind) = true;
            scantrack(ind).n = 0; % All pieces lost
            scantrack(ind).inliers = {};
        end
    end
end

%% Remove common points in different scans (corner points)
% common = cell(1,3);
% for k=1:3
%     ij = setdiff(1:3,k);
%     common{k} = intersect( new_inliers{ij(1)}, new_inliers{ij(2)} );
%     common{k} = (common{k}(:))'; % Necessary to get row vectors
% end
% common_inliers = cell2mat( common );
% for k=1:3
%     new_inliers{k} = setdiff( new_inliers{k}, common_inliers );
% end

for k=1:3
    inliers = [scantrack(k).inliers{:}]; % Stack inliers of different pieces all together
    scantrack(k).Ninliers = length(inliers);
    
    seg_pts{k} = pts(:,inliers);
end


%% Compute SVD adjustment and uncertainty propagation
Ninliers = zeros(1,3);
l = cell(1,3);
A_l = cell(1,3);
p = cell(1,3);
A_p = cell(1,3);
lin = cell(1,3);
seg = cell(1,3);
for k=1:3
    if occ(k)
        %         inliers = new_inliers{k};
        % Get complete list of inliers
        inliers = [scantrack(k).inliers{:}]; % Stack inliers of different pieces all together
        Ninliers(k) = length(inliers);
        if length(inliers) >= min_inliers
            seg_pts = pts(:,inliers);
%             [out(k), outPts{k}, lost(k)] = computeSegmentSVD(seg_pts,sigma,debug);
            [l{k}, A_l{k}, p{k}, A_p{k}, lin{k}] = computeSegmentSVD( seg_pts, sigma, debug );
            
            seg{k}(:,1) = seg_pts(:,1) - lin{k}' * makehomogeneous( seg_pts(:,1) ) * lin{k}(1:2);
            seg{k}(:,2) = seg_pts(:,end) - lin{k}' * makehomogeneous( seg_pts(:,end) ) * lin{k}(1:2);
            
            % Compute segment biggest displacement from previous iteration
%             thres_disp(k) = max( abs( scantrack(k).lin' * makehomogeneous( out(k).seg ) ) ); % To apply yet
        else
            out(k).dir = [];
            out(k).A_dir = [];
            out(k).p = [];
            out(k).A_p = [];
            out(k).line = [];
            out(k).seg = [];
            out(k).inliers = [];
            seg_pts{k} = [];
        end
    else
        out(k).dir = [];
        out(k).A_dir = [];
        out(k).p = [];
        out(k).A_p = [];
        out(k).line = [];
        out(k).seg = [];
        out(k).inliers = [];
        seg_pts{k} = [];
    end
end


%% Compute calibration data from detected segments: segment directions (l), corner points (q)
% l   = {out.dir};
% A_l = {out.A_dir};
% p   = {out.p};
% A_p = {out.A_p};
% lin = {out.line};

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

[scantrack.lin] = deal( lin{:} );
[scantrack.seg] = deal( seg{:} );

% for k=1:3
%     scantrack(k).lin = out(k).line;
%     scantrack(k).seg = out(k).seg;
% %     scantrack(k).inliers already assigned!
% 
% %     scantrack(k).n = length(scantrack(k).inliers); % Already set
% end
% Update Ninliers

% LIDAR segments should be given to this scan2co so that:
% line 1 is in plane YZ
% line 2 is in plane XZ
% line 3 is in plane XY

% Compute intersection points of three segments
% ----------------------------------------------
q = cell(1,3);

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
    plot(x,y,'.k')
    plotLIDARframe
    rgb = 'rgb';
    for k=1:3
        if ~isempty( out(k).line )
            plotHomLineWin(out(k).line,rgb(k))
            plot(out(k).seg(1,:),out(k).seg(2,:),strcat('o',rgb(k)));
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


end

