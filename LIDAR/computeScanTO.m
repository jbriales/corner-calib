function [v,A_v,p,A_p,q,A_q, lin,seg] = computeScanTO( pts, sigma, inliers, debug )
% [l,A_l,p,A_p,q,A_q, lin,seg] = computeScanCorner( scan, inliers, debug )
% Input:
%   scan - scan structure with 2D points field (xy)
%   inliers - 1x3 cell array with inlier indexes for each segment
%   debug - show debugging
% 
% Output:
%   v,A_v - line direction and 1x1 angle covariance
%   p,A_p - 2D central point and 2x2 point covariance
%   q,A_q - 2D corner point
%   lin   - 1x3 cell array with homogeneous lines through segments
%   seg   - 1x3 cell array with 2 2x1 segment limit points

% Collection of data
xy  = cell(1,3);
for k=1:3
    xy{k} = pts(:,inliers{k});
end

% Preallocation
[v,A_v,p,A_p,q,A_q,lin,seg] = deal( cell(1,3) );

%% Compute SVD adjustment and uncertainty propagation
% occ = cellfun(@(x)~isempty(x), {scantrack.lin});
occ = cellfun(@(x)numel(x)>=2, inliers);
for k=1:3
    if occ(k)
        seg_pts = xy{k};
        [v{k}, A_v{k}, p{k}, A_p{k}, lin{k}] = computeSegmentSVD( seg_pts, sigma, debug );
        
        % Compute segment end points
        seg{k}(:,1) = seg_pts(:,1) - lin{k}' * makehomogeneous( seg_pts(:,1) ) * lin{k}(1:2);
        seg{k}(:,2) = seg_pts(:,end) - lin{k}' * makehomogeneous( seg_pts(:,end) ) * lin{k}(1:2);
        
        % Compute segment biggest displacement from previous iteration
        % thres_disp(k) = max( abs( scantrack(k).lin' * makehomogeneous( out(k).seg ) ) ); % To apply yet
    else
        [v{k},A_v{k},p{k},A_p{k},lin{k},seg{k}] = deal([]);
    end
end

%% Compute calibration data from detected segments: segment directions (v), corner points (q)

% Compute homogeneous line uncertainty
A_lin = cell(1,3);
T = [ 0 -1 ; 1 0 ];
for i = 1:3
    if ~isempty( p{i} ) && ~isempty( lin{i} ) 
        % Homogeneous line computation can be done all out of SVD
%         if 0 % Perform Monte Carlo simulation
%             keyboard
%             FunMC = @(X)
%             Manifold.Dyn
%             out = Manifold.MonteCarloSim( ...,
%                 @(X)computeBackprojectedNormals( X, Rig.Camera.K ),...
%                 obj_xi, 'Ref', obj_Nbp, 'N', 1e3 );
%             keyboard
%         end
        
        % Refactor AND there was a mistake in J_l computation (depend on v)
%         J_lin_v = [0 -1; 1 0; -p{i}(2) p{i}(1)];
        J_lin_v = [ eye(2) ; -p{i}' ] * T;
%         J_lin_p = [0 0 ; 0 0; lin{i}(2) lin{i}(1)];
        J_lin_p = [ zeros(2,2) ; -(T*v{i})' ];
%         J_l = [0; 1];
        J_v_ang = T*v{i};
        A_v_ = J_v_ang * A_v{i} * J_v_ang';
        A_lin{i} = J_lin_v* A_v_ *J_lin_v' + J_lin_p* A_p{i} *J_lin_p';
        % Final result has to be rank 2
    else
        A_lin{i} = [];
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
        % Nomenclature below:
        % q  - 2D intersection point
        % qh - P2 intersection point
        % qs - spherical representation of qh
        % qr - reduced linked to S2 representation of P2 (see Forstner)

        if 1 % Without using qs
            % Compute homogeneous corner points (qh) and uncertainty
            qh = cross( lin{idx(1)}, lin{idx(2)} );
            J_qh_l1l2 = [-skew(lin{idx(2)}) skew(lin{idx(1)})];
            A_l1l2 = [A_lin{idx(1)} zeros(3,3);
                      zeros(3,3) A_lin{idx(2)}];
            J_qr_qh = null( qh' )';
            J_qh_qh  = J_qr_qh' * J_qr_qh;
            % J_qr_q is needed to drop extra rank of A_q (there is no
            % uncertainty in q direction)
            A_qh = J_qh_qh * J_qh_l1l2 * A_l1l2 * J_qh_l1l2' * J_qh_qh';
            
            % Make inhomogeneous
            q{i} = makeinhomogeneous( qh );
            J_q_qh = 1/qh(3)^2 * [ qh(3)     0 -qh(1)
                                   0     qh(3) -qh(2) ];
            A_q{i} = J_q_qh * A_qh * J_q_qh';
        else % Using qs (seems equivalent, so the other is used)
            % KEEP FOR FUTURE
            % Compute homogeneous corner points (qh) and uncertainty
            qh = cross( lin{idx(1)}, lin{idx(2)} );
            J_qh_l1l2 = [-skew(lin{idx(2)}) skew(lin{idx(1)})];
            A_l1l2 = [A_lin{idx(1)} zeros(3,3);
                      zeros(3,3) A_lin{idx(2)}]; % TOFIX
            A_qh = J_qh_l1l2 * A_l1l2 * J_qh_l1l2';
            
            % Get spherical representation and uncertainty
            qs = qh / norm(qh);
            J_qs_qh = 1/norm(qh) * (eye(3) - qs*qs');
            A_qs = J_qs_qh * A_qh * J_qs_qh';
            
            % Drop extra rank
            J_qr_qs = null( qs' )';
            J_qs_qs = J_qr_qs' * J_qr_qs;
            
            % J_qr_qs is needed to drop extra rank of A_qs (there is no
            % uncertainty in qs direction)
            A_qs = J_qs_qs * A_qs * J_qs_qs';
            
            % Make inhomogeneous
            q{i} = makeinhomogeneous( qs );
            J_q_qs = 1/qs(3)^2 * [ qs(3)     0 -qs(1)
                                   0     qs(3) -qs(2) ];
            A_q{i} = J_q_qs * A_qs * J_q_qs';
        end
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