function scantrack = manualSetScanlines( xy, mask )

plot( xy(1,:), xy(2,:), '.k' )
axis equal
plotLIDARframe

%% Manually select approximate points
% scantrack = repmat( struct('lin',[], 'seg',cell(1,0), 'inliers',cell(1,0)), 1, 3 );
scantrack = repmat( setVoidStruct( ), 1, 3 );
% lin = cell(1,3);
% seg = cell(1,3); inliers = cell(1,3);
mh = @(x) makehomogeneous(x);
label = 'XYZ';
for k=1:3
    if mask(k)
        title( sprintf('Segment on plane %s - Zoom the image - Press enter to continue',label(setdiff(1:3,k))) )
        zoom on
        pause
        title( sprintf('Segment on plane %s - Choose two points - Press enter to continue',label(setdiff(1:3,k))) )
        [px,py] = getpts;
        while mod(length(px),2) ~= 0
            warning('Nr of points has to be even, repeat selection')
            [px,py] = getpts;
        end
        user_pts = [px py]';
        
        if isempty(px)
            scantrack(k) = setVoidStruct( );
%             lin = [];
%             seg{k} = [];
        else
            nseg = numel(px)/2;
            scantrack(k).Ninliers = 0;
            for i=1:2:2*nseg
                p1 = user_pts(:,i);
                p2 = user_pts(:,i+1);
                
                % Find closest points in scan
                [p1, idx1] = minscanptdist(xy,p1);
                [p2, idx2] = minscanptdist(xy,p2);
                if idx1 > idx2 % Swap points to be in scan order
                    aux = idx2; idx2 = idx1; idx1 = aux;
                    aux = p2; p2 = p1; p1 = aux;
                    clear aux
                end
                
                lin = cross( mh(p1), mh(p2) );
                lin = lin / norm(lin(1:2));
                seg = [p1 p2];
                inliers = idx1:idx2;
                ax = axis;
                plotHomLineWin( lin, 'm' )
                plot( p1(1), p1(2), 'om' )
                plot( p2(1), p2(2), 'om' )
                axis(ax);
                % Assignment
                scantrack(k).lin = lin;
                scantrack(k).seg{end+1} = seg;
                scantrack(k).inliers{end+1} = inliers;
                scantrack(k).Ninliers = scantrack(k).Ninliers + length(inliers);
            end
            scantrack(k).n = length( scantrack(k).seg );
            scantrack(k).all_inliers = [scantrack(k).inliers{:}];
        end
    else
        scantrack(k) = setVoidStruct( );
%         scantrack(k).n = 0;
    end    
end

title( 'Final Zoom the image - Press enter to continue' )
zoom on
pause
title('Current scan segmentation')

end

function elem = setVoidStruct( )
%     elem = struct('lin',[], 'seg',cell(1,0), 'inliers',cell(1,0), 'Ninliers',0, 'n',0)
    elem.lin = [];
    elem.seg = {};
    elem.inliers = {};
    elem.Ninliers = 0;
    elem.n = 0;
    elem.all_inliers = [];
end

function [pt, idx] = minscanptdist(xy,pt)

    N = size(xy,2);
    v = xy - repmat( pt, 1, N );
    
    dist2 = sum( v.^2, 1 );
    
    [min_dist,idx] = min( dist2 ); %#ok<ASGLU>
    
    pt = xy(:,idx);
end