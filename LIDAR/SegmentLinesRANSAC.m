function S_segs = SegmentLinesRANSAC(pts,label_0,threshold,min_inliers,sigma)
% S_segs = SegmentLinesRANSAC(pts,threshold,min_inliers,sigma)
% SegmentLinesRANSAC take an input array of scan points and other
% parameters to recursively find a set of different segments in the scan,
% giving their parameters of interest

    global WITH_MONTECARLO
    if ~exist('sigma','var')
        sigma = 0.01;
    end
    if ~exist('label_0','var')
        label_0 = 0;
    end
    
    % Backup variable
%     pts_0 = pts;
    % Remove non-valid points from scan
    % TODO: Do from object construction
    pts(:, sum(pts.^2,1)==0) = [];
    Npts = size(pts,2);
    
    it = 1;    
	label = label_0;
    while it < 15 && Npts > min_inliers        
        
        [~, inliers] = ransacfitline2D(pts,threshold, 0); % Why is type single instead of double?
        
        Ninl = numel(inliers);
        if Ninl < min_inliers
            break
        end
        
        line_pts = pts(1:2,inliers);
        p = mean(line_pts, 2);
        A_p = sigma^2*eye(2) / Ninl;
        
        [~, D, V] = svd(cov(line_pts')); % Row is single observation
        v = V(:,1); % Max eigenvector, line direction
        n = V(:,2); % Min eigenvector, line normal
        l = [ n; -n'*p ];
        
        % Uncertainty propagation in direction angle
        A_pts = kron( eye(Ninl), sigma^2*eye(2) );
        J = JMan_n_p( line_pts );
        cov_angle = J * A_pts * J';
        
        %% Compute MonteCarlo simulation from X0 with estimated gaussian distribution: From oa,ob,oc and o to R
        if WITH_MONTECARLO
            Nsamples = 1e4;
            A_l = cov_angle;
            MonteCarlo_simulate( line_pts, A_pts, A_l, Nsamples ) 
        end
        
        color = 0.5 + 0.5*rand(1,3); % Brighter colours
        S_segs(it) = struct('v',   v,...
                            'A_v', cov_angle,...
                            'p',   p,...
                            'A_p', A_p,...
                            'l',   l,...
                            'eig', diag(D),...
                            'cond',D(1,1)/D(2,2),...
                            'inl', inliers,...
                            'pts', line_pts,...
                            'lab', label,...
                            'col', color); %#ok<AGROW>
                
        pts(:,inliers) = [];
        Npts = size(pts,2);
        
        it = it + 1;
        label = label + 1;
    end
end

%% Additional functions defined for uncertainty propagation and MonteCarlo simulation
function J = JMan_n_p( L0 )
% Jacobian of l tangent space angle alpha wrt set of initial points p_i
% The complete chain is:
% alpha <- l <- n <- Phi <- M <- p_i
% Input:
%   L0 is a 2xN array where each column is a point

N = size(L0,2);

c  = mean(L0,2);
L = L0 - repmat(c,1,N); % Points displaced to centroid
M = L*L'; % Scatter matrix of displaced points
[~,~,V] = svd(M);
v = V(:,1);
n = V(:,2);

% Jacobian of n wrt p_i
J_M_p = cell(1,N);
for k=1:N
    J_M_p{k} = kron(L(:,k),eye(2)) + kron(eye(2),L(:,k));
end
J_M_p = cell2mat( J_M_p );

% Implicit Function Theorem
Hess_C_nn = -n'*M'*M*n + v'*M'*M*v;
Jac_C_nM  = (M*n)' * kron( -v',eye(2) ) + ...
            (-M*v)'* kron( +n',eye(2) );
J = - Hess_C_nn \ ( Jac_C_nM * J_M_p );
end

function MonteCarlo_simulate( pts, A_pts, A_l, Nsamples ) 

X0  = [ pts(1,:) pts(2,:) ]';
F_L = @(X) adjustSegment( reshape(X,[],2)' );

fprintf('Montecarlo analysis of R (3x3) solving from Corner Observation\n')
tic, MonteCarlo.simulate( F_L, X0, 'propagation', true, 'AF', A_l, 'A0', A_pts, 'N', Nsamples, 'outSpace', 'S1' ); toc

end

function out = adjustSegment( pts )
% Function that solves direction of line (for use in MonteCarlo simulation)

N = size(pts,2);

c     = mean(pts,2);
pts_c = pts - repmat( c, 1, N );

M = pts_c * pts_c';

[~, ~, V] = svd(M);
l = V(:,1);

out = l;
end