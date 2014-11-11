function S_segs = SegmentLinesIncremental(pts,label_0,threshold,min_inliers,sigma)
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
    
    % Remove non-valid points from scan
    % TODO: Do from object construction
    pts(:, sum(pts.^2,1)==0) = [];
    
    S_segs = struct('v',   [],...
                    'A_v', [],...
                    'p',   [],...
                    'A_p', [],...
                    'l',   [],...
                    'A_l', [],...
                    'inl', [],...
                    'pts', [],...
                    'lab', [],...
                    'col', []); %#ok<AGROW>
	S_isEmpty = true;
    
    it = 1;    
	label = label_0;
    inc  = 3;
    while size(pts,2) > inc
        msd_0 = +inf;
        msd   = 0;
        msd_hist = [];
        msd_mean = +inf;
        lins = cell(1,0);
        
        tail = inc;
        % Initialization
        seg_pts = pts(:,1:tail);
        l = adjustLine( seg_pts );
        % Mean squared distance to line of support points
        d = (l' * makehomogeneous( seg_pts ))';
        msd_new = sqrt(d'*d)/numel(d);
        new_pts = [];
        while msd_new < threshold && tail+inc <= size(pts,2)
            seg_pts = [seg_pts, new_pts];
            l = adjustLine( seg_pts );
            
            new_pts = pts(:, (tail+1):(tail+inc));
            % Mean squared distance to line of new support points
            d = (l' * makehomogeneous( new_pts ))';
            msd_new = sqrt(d'*d)/numel(d);
            
            % Update and record
            tail = tail + inc;
%             msd_hist = [ msd_hist msd ];
%             lins{end+1} = l;
        end
        % Reject spurious segments
        if size(seg_pts,2) < min_inliers
            pts(:,1:tail) = [];
            continue
        end
        
        % Go back and find minimum msd from tail
%         idx = find( ([ 0 msd_hist ] - [msd_hist +inf]) > 0, 1,'last' );
%         l = lins{idx};
        % Seek inliers for this best line
        d = l' * makehomogeneous( pts );
        inliers = abs(d) < threshold;
        % Refine new line
        seg_pts = pts(:,inliers);
        [l, p, v, n] = adjustLine( seg_pts );
        
        % Uncertainty propagation
        % For central point
        Ninl = size(seg_pts,2);
        A_p = sigma^2*eye(2) / Ninl;
        % For direction
        A_pts = kron( eye(Ninl), sigma^2*eye(2) );
        J = JMan_n_p( seg_pts );
        A_angle = J * A_pts * J';
        % For homogeneous line
        % TODO: Take correlation into account
        J_l_n = [ eye(2) ; -p' ];
        J_l_p = [ zeros(2,2) ; -n' ];
        A_l = J_l_n * v * A_angle * v' * J_l_n' + ...
              J_l_p * A_p * J_l_p';

        if 0
            figure, plot(pts(1,:),pts(2,:),'.k'), hold on
            plot(seg_pts(1,:),seg_pts(2,:),'.r')
            plotHomLineWin( l, 'r' )
        end
        
        %% Compute MonteCarlo simulation from X0 with estimated gaussian distribution: From oa,ob,oc and o to R
        if WITH_MONTECARLO
            Nsamples = 1e4;
            A_l = A_angle;
            MonteCarlo_simulate( line_pts, A_pts, A_l, Nsamples ) 
        end
        
        color = 0.5 + 0.5*rand(1,3); % Brighter colours
        S_segs(it) = struct('v',   v,...
                            'A_v', A_angle,...
                            'p',   p,...
                            'A_p', A_p,...
                            'l',   l,...
                            'A_l', A_l,...
                            'inl', [],...
                            'pts', seg_pts,...
                            'lab', label,...
                            'col', color); %#ok<AGROW>
        S_isEmpty = false;
                
        % Removed used points for next segment
        pts(:,inliers) = [];
        
        label = label + 1;
        it    = it + 1;
    end
    
    if S_isEmpty
        S_segs = [];
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



function [lin, p, v, n] = adjustLine( pts )
% Function that solves direction of line (for use in MonteCarlo simulation)

N = size(pts,2);

p     = mean(pts,2);
pts_p = pts - repmat( p, 1, N );

M = pts_p * pts_p';

[~, ~, V] = svd(M);
v = V(:,1);
n = V(:,2);
lin = [ n; -n'*p ];
end