function [trel, pairs, im_points] = codeTranslation( this, segs1_all, segs2_all, matches, Rrel, K, tranThres )

    %% Inputs: Check and assign
    % Build two arrays with potentially paired segments
    if isempty(matches)
        % Skip matches argument, the two segments should be already
        % same dimension and paired by position
        if any(size(segs1_all)~=size(segs2_all))
            error('If no matches given, both segs arguments should be the same size');
        end
    else
        if size(matches,1)~=2
            warning('matches for codeP3oA should be 2xM');
            keyboard
            if size(matches,2)==2
                matches = matches';
            else
                error('Bad matches input dimension');
            end
        end
        % TEMPORAL:
        % Do not remove segments
%         segs1 = segs1_all(matches(1,:));
%         segs2 = segs2_all(matches(2,:));
        segs1 = segs1_all;
        segs2 = segs2_all;
    end
    segs = {segs1, segs2}; % To reduce code with for when possible

    % Build pairs of point data
    % TODO: Apply over all segments (out of P3oA configuration too?)
    % Create all pairs of segments (for given matches)
    N = size(matches,2);
    combs  = combnk(1:N,2)';
    Ncombs = size( combs,2 );
    
    % Construct pairs with matches indeces
    pairs = cell(2,1);
    for k=1:2
        pairs{k} = zeros(2,Ncombs);
        pairs{k}(1,:) = matches(1, combs(1,:));
        pairs{k}(2,:) = matches(1, combs(2,:));
    end
    % Store current possible triplets in CGT for later debugging
    CGT.pairs_(pairs{1}); % NOTE: Assumming pairs{1}==pairs{2}
    
    if this.TRICK_REMOVE_OUTLIERS_t % Keep only inlier configurations
        % Only work if tags or index ids are maintained!
        warning('Only GT inlier configs are being used for translation')
        pairs{1} = intersect( CGT.pairs', pairs{1}','rows' )';
        pairs{2} = intersect( CGT.pairs', pairs{2}','rows' )';
        Ncombs = size( pairs{1},2 );
    end
    
    [im_lines, im_points] = deal( cell(2,1) );
    for k=1:2
        % For each frame #1 and #2
        % Get image homogeneous lines 
        im_lines{k} = [segs{k}.l];
        % Get all combinations of these lines
        % TODO: Be careful with not-valid input segments! (Remove!?)
        % Map pair indexes to real matches indexes

        % Get image intersection point for each pair
        im_points{k} = cross( im_lines{k}(:,pairs{k}(1,:)),...
                              im_lines{k}(:,pairs{k}(2,:)), 1 );
    end
        
    % IMPORTANT: Convert points to P2 instead of image
    
    % For debug
    if 0
    figure
%     for k=1:2
%         subplot(1,2,k); hold on;
%         p{k} = makeinhomogeneous( im_points{k} );
%         plot( p{k}(1,:),p{k}(2,:),'.k');
%         segs{k}.plot('r',true);
%     end
    hold on;
    p{1} = makeinhomogeneous( im_points{1} );
    plot( p{1}(1,:),p{1}(2,:),'or');
    segs{1}.plot('r',true);
    p{2} = makeinhomogeneous( im_points{2} );
    plot( p{2}(1,:),p{2}(2,:),'*b');
    segs{2}.plot('b',true);
    end
        
    
    %% RANSAC computation
    % Apply RANSAC with epipolar constraint (from Ali)
    inliersT  = [];
    if ~exist('tranThres','var')
        % If threshold not given in input, use default value of 1[?]
        tranThres = 1; % From Elqursh (units?)
    end
    tranThresh_cur = tranThres;
    minTranInliers = min( 3, Ncombs );
    while ( length(inliersT) < minTranInliers )
        if this.WITH_DEBUG
            fprintf('\tEstimating Translation, RANSAC threshold = %f...\n',tranThresh_cur);
        end
        
        % Define parameters for RANSAC algorithm
        s = 2;
        invK = inv(K);
        fittingfn = @(pts)defineTrans(pts, Rrel', invK); % Use transpose of relative rot
        distfn    = @(t,pts,thresh)reprDist(t,pts,thresh, Rrel', invK);
        degenfn   = @isdegenerate;
        [trel, inliersT] = ransac([im_points{1};im_points{2}],...
                                  fittingfn, distfn, degenfn,...
                                  s, tranThresh_cur, false);
%         [trel,inliersT,~] = ali_t_relorient_2pts_ransac(...
%                             im_points{1},im_points{2},...
%                             K,Rrel',tranThresh_cur,... % Rrel is transposed to meet Elqursh criterion
%                             true); % True for restimate
        if this.WITH_DEBUG
            fprintf('inliers for t = %d/%d\n',length(inliersT),Ncombs);
        end
        tranThresh_cur = tranThresh_cur *2;
    end
    % Restimate with all inliers
    trel = defineTrans([im_points{1}(:,inliersT);im_points{2}(:,inliersT)],...
                        Rrel', invK);
                 
    % Given R and t find correct sign of t by testing the value of t that will
    % make the reconstructed points infront of both cameras.
    trel = ali_t_direction_chieral(Rrel',trel,...
               invK*im_points{1}(:,inliersT),...
               invK*im_points{2}(:,inliersT));
    
    % Correct trel to our criteria
    % Elqursh criterion: trel = (^c2)[(^c2)t(_c1)], translation from c2 to
    % c1, expressed in c2 coordinate system
    % P3oA criterion:    trel = (^c1)[(^c1)t(_c2)], translation from c1 to
    % c2, expressed in c1 coordinate system
    % Then, the conversion from Elqursh to ours then is:
    % (^c1)R(_c2) * [ -(^c2)[(^c2)t(_c1)] ]
    trel = -Rrel*trel;
    
    % Return inlier pairs of lines (for debug and checking)
    for k=1:2
        pairs{k} = pairs{k}(:,inliersT);
    end
    
    % Return set of inlier matched points in both images
    % A 2x1 cell array with matching between corresponding 2x1 columns (points)
    for k=1:2
        im_points{k} = im_points{k}(:,inliersT);
    end
    if 0
        % Plot inhomogeneous points
        warning('Plotting intersection points');
        k = 1; hold on, plot(im_points{k}(1,:)./im_points{k}(3,:),im_points{k}(2,:)./im_points{k}(3,:),'r*')
        k = 2; hold on, plot(im_points{k}(1,:)./im_points{k}(3,:),im_points{k}(2,:)./im_points{k}(3,:),'g*')
    end
end

function t = defineTrans(pts, R, invK)
pts1 = pts(1:3,:);
pts2 = pts(4:6,:);

% Transform
pts1 = R * invK * pts1;
pts2 = invK * pts2;

A = zeros(size(pts1,2),3);
for i=1:size(pts1,2)
    p1 = pts1(:,i);
    p2 = pts2(:,i);
    A(i,:)= [   (p2(3)*p1(2) - p2(2)*p1(3)) ...
                (p2(1)*p1(3) - p2(3)*p1(1)) ...
                (p2(2)*p1(1) - p2(1)*p1(2))];
end
[~,~,V] = svd(A);
t = V(:,end);
end

function [bestInliers, bestt] = reprDist(t,pts,thresh, R, invK)
% Square the thresh because we use the square of the distance
%         thresh = thresh^2;

%N = size(pts,2);

p1 = pts(1:3,:);
p2 = pts(4:6,:);

%p2tEp1 = zeros(1,N);

F = invK'*skew(t)*R*invK;

l1 = F'*p2;
l2 = F *p1;

% Compute the reprojection error
d = ali_dist_line2d_pnt(l1,p1) + ali_dist_line2d_pnt(l2,p2);

bestInliers = find(abs(d) < thresh);     % Indices of inlying points
bestt = t;                          % Copy F directly to bestF
end

%------------------------------------------------------------------------
% Dummy function (no degenerate sets since s=1)
function r = isdegenerate(~)
    r = false;
end