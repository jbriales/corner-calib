function [v, cov_angle, c, cov_c, lin] = computeSegmentSVD( XY, sigma, debug )
% [v,cov_angle,c,cov_c,lin] = computeSegmentSVD( XY, sigma, debug )

if ~exist('debug','var')
    debug = 0;
end

%% Compute SVD optimization and covariance propagation
pts = double( XY );
Npts = size( pts, 2 );

c = mean(pts, 2);
cov_c = sigma^2*eye(2) / Npts;

[~, ~, V] = svd(cov(pts')); % Row is single observation
v = V(:,1); % Max eigenvector, line direction
n = V(:,2); % Min eigenvector, line normal

% Uncertainty propagation in direction angle
A_pts = kron( eye(Npts), sigma^2*eye(2) );
J = J_alpha_p( pts );
cov_angle = J * A_pts * J';

% Homogeneous line
lin = [ n', -n'*c ]';

%% Compute MonteCarlo simulation from X0 with estimated gaussian distribution: From oa,ob,oc and o to R
global WITH_MONTECARLO
if WITH_MONTECARLO
    Nsamples = 1e4;
    A_l = cov_angle;
    MonteCarlo_simulate( pts, A_pts, A_l, Nsamples )
end

if debug
    fig_title = 'RANSAC Line Segmentation';
    fig = figure('units','normalized','outerposition',[0 0 1 1]); hold on, axis equal, xlabel('x'), ylabel('y'); title(fig_title)
    plot(XY(1,:),XY(2,:),'.k'); plot(0,0,'>')
%     colors =   {'g.','r.','c.','m.','y.','w.',[0.5 0.8 0.8],'k'};

    plotHomLineWin( lin, 'b' );
    axis( [ min(seg(1,:))-1, max(seg(1,:))+1, min(seg(2,:))-1, max(seg(2,:))+1 ] )
end

if debug
    keyboard
    close(fig)
end
end

%% Additional functions defined for uncertainty propagation and MonteCarlo simulation
function J = J_alpha_p( L0 )
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
l = V(:,1);
n = V(:,2);

% Jacobian of n wrt p_i
J_M_p = zeros(4,2*N);
for i=1:N
    J_M_p(:,(1:2)+(2*(i-1))) = [ 2*L(1,i)    0
                                L(2,i)  L(1,i)
                                L(2,i)  L(1,i)
                                 0     2*L(2,i) ];
end
J_Phi_M = kron( n', eye(2) );
J_n_Phi = pinv(M);
J_l_n = [0 -1; 1 0];
J_a_l = l' * [0 1 ; -1 0];
J = J_a_l * J_l_n * J_n_Phi * J_Phi_M * J_M_p;

end

function MonteCarlo_simulate( pts, A_pts, A_l, Nsamples ) 

X0  = [ pts(1,:) pts(2,:) ]';
F_L = @(X) adjustSegmentDir( reshape(X,[],2)' );

fprintf('Montecarlo analysis of R (3x3) solving from Corner Observation\n')
[manSum,manSubtraction,manMean] = MonteCarlo.spaceOps( repmat({'R'},1,size(A_pts,1)), {'S1'} );
tic, MonteCarlo.simulate( F_L, X0, A_pts, 'manSum',manSum, 'manSubtraction',manSubtraction, 'manMean',manMean,...
                          'AF', A_l, 'N', Nsamples ); toc
end

function v = adjustSegmentDir( pts )
% Function that solves direction of line (for use in MonteCarlo simulation)

N = size(pts,2);

c     = mean(pts,2);
pts_c = pts - repmat( c, 1, N );

M = pts_c * pts_c';

[~, ~, V] = svd(M);
v = V(:,1);
end