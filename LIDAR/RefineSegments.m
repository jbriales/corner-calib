function Slin = RefineSegments(seg,x,y,threshold,min_inliers,sigma,debug)
% Slin = RefineSegments(lin,x,y,threshold,min_inliers,sigma,debug)

if ~exist('debug','var')
    debug = 0;
end
if ~exist('sigma','var')
    sigma = 0.01;
end

%% Get inliers (by RANSAC) from input approximate line vector
XY = [x(:),y(:)]';
% lin = seg2line( seg );
seg0 = seg; % Store for debug

[lin, inliers] = refineLineInliers( seg, XY, threshold, min_inliers );
% Find end points in segment by projection of end inliers
seg = zeros(2,2);

seg(:,1) = XY(:,inliers(1)) - lin' * makehomogeneous( XY(:,inliers(1)) ) * lin(1:2);
seg(:,2) = XY(:,inliers(end)) - lin' * makehomogeneous( XY(:,inliers(end)) ) * lin(1:2);

%% Compute SVD optimization and covariance propagation
Ninl = length(inliers);

if Ninl >= min_inliers
    
    line_pts = double( XY(:,inliers) );
    c = mean(line_pts, 2);
    cov_c = sigma^2*eye(2) / Ninl;
    
    [~, ~, V] = svd(cov(line_pts')); % Row is single observation
    l = V(:,1); % Max eigenvector, line direction
    n = V(:,2); % Min eigenvector, line normal
    
    % Uncertainty propagation in direction angle
    A_pts = kron( eye(Ninl), sigma^2*eye(2) );
    J = J_alpha_p( line_pts );
    cov_angle = J * A_pts * J';
    
    %% Compute MonteCarlo simulation from X0 with estimated gaussian distribution: From oa,ob,oc and o to R
    global WITH_MONTECARLO
    if WITH_MONTECARLO
        Nsamples = 1e4;
        A_l = cov_angle;
        MonteCarlo_simulate( line_pts, A_pts, A_l, Nsamples )
    end
else
    warning('Min nr of inliers not found')
    l = [];
    cov_angle = [];
    c = [];
    cov_c = [];
    lin = [];
    seg = [];
    line_pts = [];
end

Slin = struct('dir', l,...
              'A_dir',cov_angle,...
              'p',   c,...
              'A_p',  cov_c,...
              'line', lin,...
              'seg', seg,...
              'pts', line_pts);

if debug
    fig_title = ['RANSAC Line Segmentation'];
    fig = figure('units','normalized','outerposition',[0 0 1 1]); hold on, axis equal, xlabel('x'), ylabel('y'); title(fig_title)
    plot(x,y,'.k'); plot(0,0,'>')
%     colors =   {'g.','r.','c.','m.','y.','w.',[0.5 0.8 0.8],'k'};

    plot(x(inliers),y(inliers),'g*')
    plotHomLineWin( lin, 'b' );
    plot(seg(1,:),seg(2,:),'or');
    plotEllipseSeg( seg0, 10*threshold, '-g' )
    plotEllipseSeg( seg,  10*threshold, '-r' )
%     plot([c(1)-2*l(1),c(1)+2*l(1)],[c(2)-2*l(2),c(2)+2*l(2)],'k-','linewidth',1)
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
% J_M_p   = chol( Mp )'; % To get equivalent jacobian from one general point
J_Phi_M = kron( n', eye(2) );
J_n_Phi = pinv(M);
J_l_n = [0 -1; 1 0];
J_a_l = l' * [0 1 ; -1 0];
J = J_a_l * J_l_n * J_n_Phi * J_Phi_M * J_M_p;

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

function [lin, inliers1] = refineLineInliers( seg, X, t, min_inliers )
% Take previous segment estimate and LIDAR scan and refine line with RANSAC

kt = 10; % Threshold scaling
if 1
    % Find points within ellipse
    [inliers1] = seg2ptDist(seg, X, kt*t);
    
else
    % Using min distance
    
    lin = seg2line( seg );
    % Threshold is bigger than in RANSAC
    [inliers1] = line2ptDist( lin, X, kt*t );
    if length(inliers1) < min_inliers
        % Get approximate set of inliers (close to previous line)
        [min_dist, pt] = minscanlindist(X,lin);
        [inliers1] = line2ptDist( lin, X, kt*min_dist);
        
        % Get refined line and inliers by RANSAC
        X1  = X( :, inliers1 );
        % lin = fitline( xy );
        [lin, ~] = ransacfitline2D(X1, t, 0);
        
        [inliers1] = line2ptDist( lin, X, kt*t );
    end
end

% Get refined line and inliers by RANSAC
X  = X( :, inliers1 );
% lin = fitline( xy );
[lin, inliers2] = ransacfitline2D(X, t, 0);
inliers1 = inliers1(inliers2); % Take good inliers of approximate set

% Drop end elements
inliers1(end-5:end) = [];
inliers1(1:5) = [];

if length(inliers1) < 5
    warning('not enough')
end

end

function [min_dist, pt] = minscanlindist(xy,lin)

    dist = abs( lin' * makehomogeneous(xy) );

    [min_dist,idx] = min( dist );
    
    pt = xy(:,idx);
end

function [inliers] = line2ptDist(l, X, t)

    mh = @(x) makehomogeneous(x);
    
    l = l / norm(l(1:2));
    
    d = abs( l' * mh(X) );
    
    inliers = find(d < t);
end

function [inliers] = seg2ptDist(seg, X, t)

    seg_d = norm( seg(:,1) - seg(:,2) );
    
    N = size( X, 2);
	seg1 = repmat( seg(:,1), 1, N );
    seg2 = repmat( seg(:,2), 1, N );
    
    % Sum of square distances
    d2 = sum( (X - seg1).^2 + (X - seg2).^2, 1);

    inliers = find( d2 <= (seg_d + 2*t)^2 );
end

function plotEllipseSeg( seg, t, format )
d = norm( seg(:,1) - seg(:,2) );

b = d/2 + t;
a = sqrt( b^2 - (d/2)^2 );

th = 0:0.1:2*pi+0.1;
ct = cos( th );
st = sin( th );

pts = [ a * ct ; b * st ];
c = mean( seg, 2 );
alpha = atan2( seg(2,2) - seg(2,1), seg(1,2) - seg(1,1) ) + pi/2;
R = [ cos(alpha) -sin(alpha)
      sin(alpha)  cos(alpha) ];
N = length( th );
pts = repmat(c, 1, N) + R * pts;
plot( pts(1,:), pts(2,:), format )

end

function lin = seg2line( seg )
lin = cross( makehomogeneous( seg(:,1) ), makehomogeneous( seg(:,2) ) );
lin = lin / norm(lin(1:2));
end