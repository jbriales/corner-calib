function [ R, cov, cov_eps, err, errNorm, W ] = optimRotation( rot_input, R, weighted, all_label, gt )
% R = optimRotation( rot_input, R, weighted, gt )
% 
% rot_input is a array of structs where each element of the array
% corresponds to an single observation of a corner and each struct
% contains:
%   N - normal vectors to world planes
%   A_N - uncertainty of N
%   A_eps - uncertainty of eps linked to every N triple (for Monte Carlo)
%   l - direction vector of segments in LIDAR SR
%   A_l - uncertainty of l
% This format allows to process the set of observations as blocks (each
% observation could have 3 or less pieces of data)

global WITH_MONTECARLO

if ~exist('gt','var')
    gt = [];
end

n_cam = [rot_input.N];
line_LRF = [rot_input.l];
if ~weighted
    Lev_Fun = @(R) Fun( n_cam, line_LRF, R );
    [ R, err, errNorm, W ] = LM_Man_optim(Lev_Fun,R,'space','SO(3)','debug',2, 'maxIters', 200);
    
    % Compute R covariance
    Jac = optimDer( rot_input, R );
    % Create mapping and cov
    N = [rot_input.N];
    A_N = blkdiag(rot_input.A_N);
    L = [rot_input.l];
    A_L = blkdiag(rot_input.A_l);
    A = blkdiag( A_N, A_L );
    
    % Permutation matrix for n(3x1) and l change in Lie algebra (1x1)
    T = zeros(4*size(N,2));
    for i=1:size(N,2)
        for k=1:3
            abk = (i-1)*4 + k;
            abv = (i-1)*3 + k;
            T(abk,abv) = 1;
        end
        for k=1:1
            abk = 3 + (i-1)*4 + k;
            abv = numel(N) + (i-1)*1 + k;
            T(abk,abv) = 1;
        end
    end
    A_eps = Jac * T * A * ( Jac * T )';
    J_R_eps = MonteCarlo.manDiffRotExp(R);
	A_R = J_R_eps * A_eps * J_R_eps';
    cov_eps = A_eps;
    cov = reshape( diag(A_R), 3, 3 );
    
    if WITH_MONTECARLO
        MonteCarlo_simulate( );
    end
else
    Lev_Fun = @(R) FunW( rot_input, R );
    [ R, err, errNorm, W ] = LM_Man_optim(Lev_Fun,R,'space','SO(3)','weighted',true,'debug',2, 'maxIters', 200);   
    
    cov_eps = [];
    cov = []; % TODO: Compute covariance
end

    function MonteCarlo_simulate( )
        X0 = [ n_cam(:); line_LRF(:) ];
        function R = F_optim(X0)
            n = length( X0 ) / 5;
            N = reshape(X0(1:n*3),3,n);
            L = reshape(X0(n*3+1:end),2,n);
            Lev_Fun = @(R) Fun( N, L, R );
            [ R, ~,~,~ ] = LM_Man_optim(Lev_Fun,R,'space','SO(3)','debug',2, 'maxIters', 200);
        end

        A_N_man = blkdiag(rot_input.A_eps);
        % TODO: Make new manifold rule for normals pair        
        % Create list of inputs:
        Cin = cell(1,length(rot_input)+size(line_LRF,2)); % TODO: Correct
        for k=1:length(rot_input)
            switch size(rot_input(k).N,2)
                case 1
                    man = struct('topo','R','dim',3);
                case 2
                    man = struct('topo','SO*','dim',3);
                case 3
                    man = struct('topo','SO','dim',3);
            end
            Cin{k} = man;
        end
        [Cin{length(rot_input)+1:end}] = deal( struct('topo','S','dim',1) );
        Cout = {struct('topo','SO','dim',3)};
        
        fprintf('Montecarlo analysis of R covariance from N and L\n')
        [manSum,manSubtraction,manMean] = MonteCarlo.spaceOps( Cin, Cout );
        tic, [mu,A] = MonteCarlo.simulate( @F_optim, X0, blkdiag(A_N_man,A_L), 'manSum',manSum, 'manSubtraction',manSubtraction, 'manMean',manMean,...
            'AF', A_eps, 'N', 1e4 ); toc
        keyboard
    end

% Plotting results and weights
figure, hold on
bar_width = 0.25;
title(sprintf('Optimization with WEIGHTED=%d',weighted))
all_label( all_label==0 ) = [];
rgb = 'rgb';
subplot(411), hold on
% for i=1:length(all_label)
%     bar(i,err(i), rgb(all_label(i)))
% end
for i=1:3
    bar( find(all_label==i), err(all_label==i), bar_width, rgb(i) )
end

subplot(412), hold on
errW = W*err;
for i=1:3
    bar( find(all_label==i), errW(all_label==i), bar_width, rgb(i) )
end

subplot(413), hold on
% for i=1:length(all_label)
%     bar(i,W(i,i), rgb(all_label(i)))
% end
for i=1:3
    v = diag(W);
    bar( find(all_label==i), v(all_label==i), bar_width, rgb(i) )
end
subplot(414), hold on
axis ij, colormap(gray)
handle = bar3(abs(W)); % abs to set 0 as black
colorbar
for k=1:length(handle)
    zdata = get(handle(k),'ZData');
    set(handle(k),'CData',zdata,'FaceColor','interp');
end



end

function [residual, J] = Fun( n_cam, line_LRF, R )

% n_cam = [rot_input.N];
% line_LRF = [rot_input.l];
residual = dot( n_cam, R(1:3,1:2) * line_LRF, 1 )';
J = cross( R(1:3,1:2) * line_LRF, n_cam, 1 )';

end

function [residual, J, W] = FunW( rot_input, R )
global WITH_MONTECARLO

n_cam = [rot_input.N];
line_LRF = [rot_input.l];
residual = dot( n_cam, R(1:3,1:2) * line_LRF, 1 )';
J = cross( R(1:3,1:2) * line_LRF, n_cam, 1 )';

% Weights
N_obs = length(rot_input);
W = cell(1, N_obs);

ort = [ 0 -1
        1  0 ];
    
R12 = R(1:3,1:2); % Only 2 columns are used

for i=1:N_obs
    N = rot_input(i).N;
    A_N = rot_input(i).A_N;
    L = rot_input(i).l;
    A_L = rot_input(i).A_l;
        
    s  = size( N, 2 );
    RL = mat2cell( R12*L, 3, ones(1,s) );
    JN = blkdiag( RL{:} )';
    Ja = diag( dot( N, R12*ort*L ) );
    A_k = JN * A_N * JN' + Ja * A_L * Ja';

    W{i} = pinv(A_k);
    if WITH_MONTECARLO
        if s==3
            MonteCarlo_simulate( rot_input(i).A_eps )
        end
    end
end
% Build complete matrix from diagonal blocks
W = blkdiag(W{:});
% % Temporal debug:
% W = diag(diag(W));

    function MonteCarlo_simulate( A_eps )
        
        X0 = [ N(:); L(:) ];
        function e = F_e(X0)
            N = reshape(X0(1:3*3),3,3);
            L = reshape(X0(3*3+1:end),2,3);
            e = dot( N, R(1:3,1:2)*L, 1 )';
        end

        Cin = cell(1,1+size(L,2)); % TODO: Correct
        switch size(N,2)
            case 1
                man = struct('topo','R','dim',3);
            case 2
                man = struct('topo','SO*','dim',3);
            case 3
                man = struct('topo','SO','dim',3);
        end
        Cin{1} = man;
        [Cin{2:end}] = deal( struct('topo','S','dim',1) );
        Cout = {struct('topo','R','dim',3)};
        
        fprintf('Montecarlo analysis of R covariance from N and L\n')
        [manSum,manSubtraction,manMean] = MonteCarlo.spaceOps( Cin, Cout );
        tic, [mu,A] = MonteCarlo.simulate( @F_e, X0, blkdiag(A_eps,A_L), 'manSum',manSum, 'manSubtraction',manSubtraction, 'manMean',manMean,...
            'AF', A_k, 'N', 1e5 ); toc
        keyboard
    end

end


function J = optimDer( rot_input, R )
% Compute Jacobian of R wrt input parameters n and l by implicit function
% theorem
n_cam = [rot_input.N];
line_LRF = [rot_input.l];

N = size(n_cam,2); % Nr of pairs
% H = zeros(3,3);
% D = zeros(3,5*N);
blockH = cell(1,N);
blockD = cell(1,N);
% Compute Hessian and 2nd der
for i=1:N % For each pair
    ni  = n_cam(:,i); % Normal vector seen from camera (3x1)
    lsi = line_LRF(:,i); % Line vector seen from scanner (2x1)
    lci = R(:,1:2) * lsi; % Line vector seen from camera (3x1)
    
    % Phi is dC/deps = -Sum( cross(R*lsi, ni*(R*lsi)'*ni) )
    d_phi_ni  = 2 * skew(lci) * (ni*lci' + lci'*ni*eye(3));
    d_phi_lci = 2 * skew(lci) * (ni*ni') - 2 * skew(ni * lci' * ni);
    d_lci_lsi = R(:,1:2);
    d_lci_eps = - skew( R(:,1:2)*lsi );
    
    % Derivation wrt l Lie algebra instead of complete space
    d_lsi_alpha = [ -lsi(2), +lsi(1) ]';
    
%     H = H + d_phi_lci * d_phi_eps;
    blockH{i} = d_phi_lci * d_lci_eps;
%     blockD{i} = [ d_phi_ni , d_phi_lci * d_lci_lsi ];
    blockD{i} = [ d_phi_ni , d_phi_lci * d_lci_lsi * d_lsi_alpha ];
end
H = sum( reshape(cell2mat( blockH ),3,3,[]), 3 );
D = cell2mat( blockD );

J = - H \ D; % Jacobian of rotation in Lie algebra wrt all inputs ni and li
end