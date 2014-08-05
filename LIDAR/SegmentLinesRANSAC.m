function Slin = SegmentLinesRANSAC(x,y,threshold,min_inliers,sigma,debug)

    if ~exist('debug','var')
        debug = 0;
    end
    if ~exist('sigma','var')
        sigma = 0.01;
    end
    
    xyz = [x;y;zeros(1,size(x,2))];
    if debug
        fig_title = ['RANSAC Line Segmentation'];
        fig = figure('units','normalized','outerposition',[0 0 1 1]), hold on, axis equal, xlabel('x'), ylabel('y'), title(fig_title)
        plot(xyz(1,:),xyz(2,:),'.'); plot(0,0,'>')      
        colors =   {'g.','r.','c.','m.','y.','w.',[0.5 0.8 0.8],'k'};        
    end
    
    it = 1;
    remaining_pts = size(x,2);
    Slin = struct('dir',[],...
                  'A_dir',  [],...
                  'p',   [],...
                  'A_p',  [] );
    while it < 5 && remaining_pts > min_inliers        
        
        [~, ~, inliers] = ransacfitline(xyz,threshold, 0); % Why is type single instead of double?
        
        if size(inliers,1) < min_inliers
            break
        end
        Ninl = size(inliers,1);
        
        line_pts = double( xyz(1:2,inliers) );
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
        
        Slin(it) = struct('dir', l,...
                          'A_dir',cov_angle,...
                          'p',   c,...
                          'A_p',  cov_c );

    
        if debug
            %plot(V(1,1:2),V(2,1:2),'r*')
            plot(xyz(1,inliers),xyz(2,inliers),colors{it})
            plot([c(1)-2*l(1),c(1)+2*l(1)],[c(2)-2*l(2),c(2)+2*l(2)],'k-','linewidth',1)
        end  
        
        remaining_pts = remaining_pts - size(inliers,1);
        xyz = xyz(:,setdiff(1:size(xyz,2),inliers)');
        
        it = it + 1;
    end
    
    % Arrange segments wrt Y-coordinate (negative to positive -> right to
    % left)
    p = [Slin(:).p];
    [~,ind_l] = sort(p(2,:));
    Slin = Slin(ind_l);
    
    if debug
        pause;
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