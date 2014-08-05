function [scantrack, seg_pts, lost] = updateSegments( scan, scantrack, debug )
% [scantrack, seg_pts, lost] = updateSegments( scan, scantrack, debug )
% Input:
%   scan - scan structure with field xy
%   scantrack - scan tracking structure
% Output:
%   scantrack - updated structure for tracking
%   seg_pts - 1x3 cell with established inliers for each segment
%   lost - 1x3 logical array with just lost segments

if ~exist('debug','var')
    debug = false;
end

% Control variables
threshold = 0.02;
min_inliers  = 30;
min_inliers2 = 5; % For secundary pieces of line
% sigma = 0.01; % SD of LIDAR sensor: Hokuyo UTM-30LX: 0.010 <10m and 0.030 <30m
% debug = 0;

% Input data
pts = scan.xy;

seg_pts = cell(1,3);
lost = false(1,3);
lin = cell(1,3);

%% Update segments and inliers
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
            if scantrack(k).n == 0
                lost(k) = true;
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
    scantrack(k).all_inliers = inliers;
    scantrack(k).Ninliers = length(inliers);
    
    seg_pts{k} = pts(:,inliers);
end