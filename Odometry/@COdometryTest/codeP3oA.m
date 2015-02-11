function [Rrel, Vs, triplets] = codeP3oA( this, segs1_all, segs2_all, matches, K, rotThres )

% Currently segs1 and segs2 have same dimension and are paired, later this
% could be generalized with matches array
    
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
    
    %% Compute rotation
    % Create all possible combinations of three different segments in im1
    N = size(matches,2);
    triplets  = combnk(1:N,3)';
    Ntriplets = size( triplets,2 );
    % TEMPORAL: Map simple indexes to original all_segments indexes
    % Elqursh code: matches is Mx2, seg_pairs Nx2
%     seg_pairs1 = matches(seg_pairs);
%     m2 = zeros(1,size(matches,1));
%     m2(matches(:,1)) = matches(:,2);
%     seg_pairs2 = m2(seg_pairs1);
    for row=1:3
        triplets(row,:) = matches(1,triplets(row,:));
    end
    
    % Store current possible triplets in CGT for later debugging
    CGT.triplets_(triplets);
       
    % Filter impossible configurations
    mask_feasible = true(1,Ntriplets); % True for feasible configurations
    for i=1:Ntriplets
        if 0
            % Case of extension method
            % Check signature of associated quadratic form Q
            s = zeros(1,2);
            % Check interpretation angles condition: True if violated
            bp_ang_cond = false(1,2);
            for k=1:2 % For each image triplet
                Nbp  = snormalize( K' * [segs{k}( triplets(:,i) ).l] );
                s(k) = CP3oASolver.signature( Nbp );

                if CP3oASolver.isMeeting( Nbp )
                    if CP3oASolver.cos_prod( Nbp ) <= 0
                        % Equivalent to?: Check minimum angles constraint (sum > 180)
                        bp_ang_cond(k) = true;
                    end
                end
            end

            if any( s>0 ) || any( bp_ang_cond )
                mask_feasible(i) = false; % Write false if any condition is not fulfilled
            end
        else
            % Case of general method
            A = zeros(1,2);
            for k=1:2
                Nbp  = snormalize( K' * [segs{k}( triplets(:,i) ).l] );
                A(k) = CP3oASolver.discriminant( Nbp ); % Check matrix discriminant
            end
%             if any( A <= 0 )
            if any( A < CP3oASolver.thresDiscr ) % HEURISTIC thresDiscr is a constant stored in class
                mask_feasible(i) = false;
            end
        end
    end
    
    triplets = triplets(:,mask_feasible); % Keep only good ones
    Ntriplets = size( triplets,2 );
    
    if 0 % Remove outlier configurations
        % Only work if tags or index ids are maintained!
        disp('Outliers are being cropped using GT in codeP3oA for rotation')
        out_triplets = setdiff(triplets',CGT.triplets','rows')';
        triplets = setdiff( triplets', out_triplets','rows')';
        Ntriplets = size( triplets,2 );
    end
    if this.TRICK_REMOVE_OUTLIERS_R % Keep only inlier configurations
        % Only work if tags or index ids are maintained!
        warning('Only GT inlier configs are being used for P3oA Rotation')
        triplets = intersect( CGT.triplets', triplets','rows' )';
        Ntriplets = size( triplets,2 );
    end
    
    % Find 3D lines directions for each triplet in each image
    V   = cell(1,Ntriplets);
    A_V = cell(1,Ntriplets);
    Vlin = zeros( 4*3*3, Ntriplets );
    
    % If there are no triplets left after filtering, return NaN
    if Ntriplets == 0
        Rrel = NaN(3);
        Vs = {[],[]};
        triplets = [];
        return
    end
    
    % Build solver (methods factory, only keep K)
    P3oA_solver = CP3oASolver( K );
    % TODO: Improve speed by reducing idle operations
    for i=1:Ntriplets % For each triplet
        % Solve orthogonal directions for each triplet in both images and
        % make signs of directions consistent
        % NOTE: V_{i}{j} - i is the Necker configuration, j the frame number
        V_ =   { cell(2,1); cell(2,1) }; % Temporal directions cell
        A_V_ = { cell(2,1); cell(2,1) }; % Temporal covariance cell
        for k=1:2 % For each image
            Nbp  = snormalize( K' * [segs{k}( triplets(:,i) ).l] );
                        
%             out = P3oA_solver.solve( Nbp );
            out = P3oA_solver.generalSolve( Nbp );
            V_{1}{k} = out{1};
            V_{2}{k} = out{2};
            
            % Compute covariance for solutions
            % EXPERIMENTAL
            if this.WITH_COVARIANCE
                A_L = { segs{k}( triplets(:,i) ).A_l }; % Collect list output of array CSegment2D get method
                A_N = P3oA_solver.covProp_N_L( A_L, segs{k}( triplets(:,i) ) );
                A_V_{1}{k} = P3oA_solver.covProp_eps_N( A_N, Nbp, V_{1}{k} );
                A_V_{2}{k} = P3oA_solver.covProp_eps_N( A_N, Nbp, V_{2}{k} );
            end
        end
        
        % In general case two possible solutions (with undetermined sign)
        % are obtained for each triplet image
        if 0
            % We seek the good pairs. Heuris: since signs are not
            % determined, compare Frobenious norm of abs(V)
            froDist = zeros(2,2);
            % froDist: row ii has absolute Fro distance of ii-th V1 sol to both
            % possible V2 solutions
            index2 = zeros(2,1);
            % index2: ii-th element contains index of good solution in V{2}{jj}
            % for V{1}{ii}
            for ii=1:2
                for jj=1:2
                    froDist(ii,jj) = norm( abs(V_{ii}{1}) - abs(V_{jj}{2}), 'fro' );
                end
                [~,index2(ii)] = min(froDist(ii,:),[],2);
            end
        else
            % We seek the good pairs. Heuris: since signs are not
            % determined, compare abs of dot products
            % between corresponding directions (should be ~1) [dot-test]
            dotSum = zeros(2,2);
            % dotSum: row ii has sum of absolute dot product
            % of ii-th V1 conf to both possible V2 Necker configurations
            index2 = zeros(2,1);
            % index2: ii-th element contains index of good configuration in V{jj}{2}
            % for V{ii}{1}
            for ii=1:2
                for jj=1:2
                    dotSum(ii,jj) = sum( abs( dot( V_{ii}{1}, V_{jj}{2} ) ) );
                end
            end
            
            % Take the two most likely pairs of configurations
            % based on the abs-dot test
            [Y,I] = sort(dotSum(:),1,'descend');
            % If max Y value is significantly lower than 3,
            % it could be because:
            % - The true rotation angle is not low
            % - This segments pair is an outlier
% % % %           TO REMOVE
% % %             if Y(1) < 2
% % %                 warning('Very low dot matching value, check');
% % %                 keyboard
% % %             end
            [ii,jj] = ind2sub( [2 2], I(1:2) );
            % ii and jj contain, respectively, the conf index to take
            % from frames #1 and #2 configurations
            % in the corresponding new most-likely configuration
            temp_ = V_;
            for k=1:2
                V_{k}{1} = temp_{ii(k)}{1};
                V_{k}{2} = temp_{jj(k)}{2};
            end
            
            % Correct sign of both direction matrices
            for n=1:2 % For each Necker configuration pair
                V_{n}{2} = V_{n}{2} * diag( sign( dot( V_{n}{1}, V_{n}{2} ) ) );
            end
        end
        
        for k=1:2
            % TODO: Make some kind of checking in values? abs(...)~1
            % Here we assume that the Rotation is at most 90 degrees
            % Choose column signs for V{...} such that it gives
            % a positive dot product between corresponding directions
            S = diag( sign( dot( V_{k}{1}, V_{k}{2}, 1 ) ) );
            V_{k}{2} = V_{k}{2} * S;
        end
                
        % Store directions as triplet directions
        V{i} = V_;
        % Store each double-double triplet (4 sets of 3D vectors)
        % linearized for RANSAC algorithm
        % For each column:
        % 1-9:   Conf #1 in image #1
        % 10-18: Conf #1 in image #2
        % 19-27: Conf #2 in image #1
        % 28-36: Conf #2 in image #2
        Vlin(:,i) = [ V_{1}{1}(:);
                      V_{1}{2}(:);
                      V_{2}{1}(:);
                      V_{2}{2}(:) ];
                  
        % Compute and store covariance or uncertainty of each result
        if this.WITH_COVARIANCE
%             A_V{i} = P3oA_solver.computeWeight( A_N, Nbp, V_tri );
        end
    end

    % Compute all rotations and analyze statistically or with RANSAC
    all_R = cell(2,1);
    for n=1:2
        % REMEMBER:
        % For each column:
        % 1-9:   Conf #1 in image #1
        % 10-18: Conf #1 in image #2
        % 19-27: Conf #2 in image #1
        % 28-36: Conf #2 in image #2
        % Vectorized directions matrix in images #1 and #2, config #n
        V1 = Vlin((n-1)*18+ (1:9),  :);
        V2 = Vlin((n-1)*18+ (10:18),:);
        
        % Compute all rotations R{m} = V1*V2'
        % NOTE: Using manopt toolbox for multiarray operations
        all_R{n}  = multiprod( reshape(V1,3,3,[]),...
                            multitransp(reshape(V2,3,3,[])) );
    end
    Rlin = [ reshape(all_R{1},9,[])
             reshape(all_R{2},9,[]) ]; % Linealized results for RANSAC
	RRlin = reshape(Rlin,9,[]);
    
    R_GT = CGT.R;
    all_R_err_GT = [ 2 * asind( 1/(2*sqrt(2)) * sqrt(sum(...
                    (reshape(all_R{1},9,[]) - repmat(R_GT(:),1,size(all_R{1},3))).^2, 1 ) ));...
                  2 * asind( 1/(2*sqrt(2)) * sqrt(sum(...
                    (reshape(all_R{2},9,[]) - repmat(R_GT(:),1,size(all_R{2},3))).^2, 1 ) )) ];
%     figure,hist(all_R_err', [0.1 0.5 1 2 5 10]);
%     figure,hist(all_R_err_GT')
%     figure,hist(all_R_err',[0:0.1:2.5])
    
    if this.WITH_KDE
    % Try KDE for estimating rotation components and find matrix
    % close to solution filtering outliers
%     tic
    figure
    for m = 1:9
        if 0
        data = RRlin(m,:);
        [bw,pdf,xmesh,cdf] = kde(data,2^14,min(data)-0.1,max(data)+0.1);
%         [bw,pdf,xmesh,cdf] = kde(data,2^8,min(data),max(data)); %#ok<ASGLU>
        % Find max mode value
        [~,I] = max(pdf);
        R_(m) = xmesh(I);
        subplot(3,3,m);
        kde(data,2^8,min(data),max(data));
        else
%         [pdf,xi,bw] = ksdensity( RRlin(m,:) );
        % Compute heuristic bw from median distance between samples
        % in order to get continuity in most dense area only
        % TODO: Improve speed of this code
%         median_dist = median( diff(sort(RRlin(1,:))) );
%         bw = median_dist * 1000;
        bw = 0.01; % Use proportional to optimal bw? Must be as low as possible
        % to compute mode with precision
        [pdf,xi,bw] = ksdensity( RRlin(m,:), 'bandwidth', bw, 'support', [-1 1] );
        [~,I] = max(pdf);
        R_(m) = xi(I);
        subplot(3,3,m);
        plot(xi,pdf);
        end
    end
    R_ = reshape(R_,3,3);
    [Ur,~,Vr] = svd( R_ ); R = Ur*Vr'; clear Ur Vr;
%     toc
    
    all_R_err = [ 2 * asind( 1/(2*sqrt(2)) * sqrt(sum(...
                    (reshape(all_R{1},9,[]) - repmat(R(:),1,size(all_R{1},3))).^2, 1 ) ));...
                  2 * asind( 1/(2*sqrt(2)) * sqrt(sum(...
                    (reshape(all_R{2},9,[]) - repmat(R(:),1,size(all_R{2},3))).^2, 1 ) )) ];
    
    % Use error wrt estimated R to choose inliers
    % NOTE: kde (Kernel Density Estimation) could be used to automate threshold selection, but more
    % insight would be necessary in kde first
%     data = all_R_err(:);
%     [bw,pdf,xmesh,cdf] = kde(data,2^8,min(data),max(data)); %#ok<ASGLU>
%     kde(data,2^14,0,max(data));
    % Estimate rotation inliers
    % The err metric is the distance angle between rotations
    if ~exist('rotThres','var')
        % If threshold not given in input, use default value of 3[deg]
        rotThres = 1.5; % In degrees
    end
    % Choose minor error for both available Necker configurations
    [angDist, iNecker] = min( all_R_err,[],1 );
    inliers = find( angDist < rotThres );
    
    % If no inliers left skip next steps and return NaN
    if numel(inliers) == 0
        warning('No inliers found');
%         keyboard
        Rrel = NaN(3);
        Vs = {[],[]};
        triplets = NaN(3,1);
        return
    end
    
    Vinl = [V{:}];
    NV = size(Vinl,2);
    linidx = sub2ind([2,NV], iNecker(inliers), inliers);
    Vinl = [Vinl{linidx}];
    V1 = [Vinl{1,:}];
    V2 = [Vinl{2,:}];
    % Solution to ||D1-R'*D2|| is SVD of D2*D1'
    % so solution to ||D1-R *D2|| is SVD of D1*D2'
    [Ur,~,Vr] = svd(V1*V2');
    if (det(Ur)*det(Vr)>=0)
        Sr = eye(3);
    else
        Sr = diag([1 1 -1]);
    end
    Rrel = Ur*Sr*Vr';
    
    if nargout > 1
        % If output arrays of directions are demanded
        Vs = { V1; V2 };
    end
    % Keep only inlier triplets
    triplets = triplets( :, inliers );
    
    % Exit here for current implementation (not using RANSAC)
    % The code below is kept for later 
    return
                    
    
    else % WITH RANSAC instead of WITH_KDE
    % Estimate rotation inliers
    % For now, the err metric is the sum of angles between all 3 directions
    if ~exist('rotThres','var')
        % If threshold not given in input, use default value of 3[deg]
        rotThres = 1.5; % In degrees 
    end
    rotThres_cur = rotThres;
    minRotInliers = min( 3, Ntriplets );
    inliers = struct('idxs', [],...
                     'conf', [],...
                     'err',  []);
%     tic
    counter = 10;
    while numel(inliers.idxs) < minRotInliers && counter > 0
        if this.WITH_DEBUG
            fprintf('\tEstimating Rotation , RANSAC threshold = %f ...\n',rotThres_cur);
        end
        
        % Define parameters for RANSAC algorithm
        s = 1;
        fittingfn = @defineRot;
        distfn    = @angularDist;
        degenfn   = @isdegenerate;
        [Rrel, inliers] = ransac(Rlin, fittingfn, distfn, degenfn,...
                              s, rotThres_cur, this.WITH_DEBUG);
        % In case of failing RANSAC:
        if isempty(inliers)
            inliers.idxs = [];
        end
        if this.WITH_DEBUG
            fprintf('Inliers for R = %d/%d\n', numel(inliers.idxs),Ntriplets);
        end
        rotThres_cur = rotThres_cur *2;    
        counter = counter - 1;
    end
    % Undo last step
    rotThres_cur = rotThres_cur /2;
%     toc

    if isempty(inliers.idxs)
        Rrel = NaN(3,3);
        Vs = { NaN; NaN };
        return
    end

    % Restimate R from all inliers
    % using Procrustes (minimize |v1 - Rrel*v2|)
    % Do in a while loop so that the process of inliers selection can be
    % done more robustly using better estimates of R
    if (numel(inliers.idxs) > 1)
        counter = 1;
        while counter <= 10
            % Choose inliers
            Vinl = [V{inliers.idxs}];

            % Choose good configurations
            mask_of_Necker_conf = inliers.conf( inliers.idxs );
            linear_idxs = sub2ind( [2 numel(inliers.idxs)],...
                mask_of_Necker_conf, 1:numel(inliers.idxs) );
            Vinl = [Vinl{linear_idxs}];
            V1 = [Vinl{1,:}];
            V2 = [Vinl{2,:}];
            % Solution to ||D1-R'*D2|| is SVD of D2*D1'
            % so solution to ||D1-R *D2|| is SVD of D1*D2'
            [Ur,~,Vr] = svd(V1*V2');
            if (det(Ur)*det(Vr)>=0)
                Sr = eye(3);
            else
                Sr = diag([1 1 -1]);
            end
            Rrel = Ur*Sr*Vr';

            % Filter inliers again with new estimate
            % NOTE: Using original rotThres
            [inliers_test, ~] = angularDist({Rrel}, Rlin, rotThres_cur);

            if ~isequal(inliers_test.idxs, inliers.idxs)
                % If inliers change, repeat refinement process with new inliers
                inliers.idxs = inliers_test.idxs;
                counter = counter + 1;
            else
                % If inliers are the same, exit refinement loop
                break
            end
        end
    end
    % Finally copy into main variable
    inliers.idxs = inliers_test.idxs;
    if this.WITH_DEBUG
        fprintf('%d refinement loops have been done\n',counter);
        fprintf('Inliers for R = %d/%d\n', numel(inliers_test.idxs),Ntriplets);
    end
    
    if nargout > 1
        % If output arrays of directions are demanded
        Vs = { V1; V2 };
    end
    
    % Give commited error for each inlier
    % NOTE: Using original rotThres
    [inliers_test, R_test] = angularDist({Rrel}, Vlin(:,inliers.idxs), rotThres_cur); %#ok<NASGU,ASGLU>
    % Debug
%     [inliers_test, R_test] = angularDist({CGT.R}, Vlin(:,inliers.idxs), rotThres);
%     [inliers_test, R_test] = angularDist({Rrel},  Vlin, rotThres)
%     [inliers_test, R_test] = angularDist({CGT.R}, Vlin, rotThres)

    triplets = triplets( :, inliers.idxs );
    end
end

function cell_R = computeRot( cell_V )
    % Take the two possible solutions
    cell_R = { cell_V{1}{1} / cell_V{1}{2};...
               cell_V{2}{1} / cell_V{2}{2} };
end

%------------------------------------------------------------------------
% Function to define one or two possible rotations given a linearized
% a set of 3D directions in both images.
function R = defineRot(X)
    % Reshape received column data
    R = cell(1,2);
    R{1} = reshape(X(1:9),3,3);
    R{2} = reshape(X(10:18),3,3);
end
    
%------------------------------------------------------------------------
% Function to calculate angular distances between rotation candidates.
% Ideally all inlier rotations should have the same expression.
function [inliers, R] = angularDist(R, X, t)

    cell_inliers = cell(numel(R),1); % Cell for inlier indexes
    cell_d = cell(numel(R),1); % Cell for angular errors per triplet
    cell_I = cell(numel(R),1); % Cell for good configuration index (solve Necker cube)  

    NX = size(X,2);
    for n=1:numel(R)
        % For each input Rotation potential candidate (usually two, but
        % could be one)
        RR = repmat( R{n}(:), 1, NX );
        angDist_ = 2 * asind( 1/(2*sqrt(2)) * ...
                        [ sqrt(sum((X(1:9,:)-RR).^2,1))
                          sqrt(sum((X(10:18,:)-RR).^2,1)) ] );
        % Keep only lower error of both possible configurations
        [cell_d{n},cell_I{n}] = min( angDist_,[],1 );
        cell_inliers{n} = find(cell_d{n} < t);
    end
    
    % Choose best solution according to number of inliers and global error
    % HEURISTIC:
    % 1. Weight number of inliers with committed error to compare
    % both options
    % 2. Reject solutions with less than a minimum number of inliers,
    % to avoid kind of over-adjustment.
    if numel(R)>1 % Case two solutions are being evaluated
        score = zeros(1,2);
        for n=1:2
            if numel(cell_inliers{n}) <= 2 % HARD-CODED
                cell_inliers{n} = [];
            end
            % Score: kind of mean error?
%             score(n) = sum( cell_d{n}( cell_inliers{n} ) ) / ...
%                        numel(cell_inliers{n})^2; % Squared to weight numel importance
%             score_(n) = median( cell_d{n}( cell_inliers{n} ) );
            score(n) = median( cell_d{n}( cell_inliers{n} ) ) /...
                       numel(cell_inliers{n});
        end    
        [~,best_idx] = min( score );
    else
        best_idx = 1; % The only case available
    end
    % Take best candidate
    R = R{best_idx};
    % inliers = cell_inliers{best_idx};
    inliers = struct('idxs', cell_inliers{best_idx},...
                     'conf', cell_I{best_idx},...
                     'err',  cell_d{best_idx});
end
    
%------------------------------------------------------------------------
% Dummy function (no degenerate sets since s=1)
function r = isdegenerate(~)
    r = false;
end