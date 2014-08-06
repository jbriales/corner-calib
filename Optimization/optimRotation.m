function [ R, err, errNorm, W ] = optimRotation( rot_input, R0, weighted, all_label, gt )
% R = optimRotation( rot_input, R0, weighted, gt )
% 
% rot_input is a array of structs where each element of the array
% corresponds to an single observation of a corner and each struct
% contains:
%   N - normal vectors to world planes
%   A_N - uncertainty of N
%   l - direction vector of segments in LIDAR SR
%   A_l - uncertainty of l
% This format allows to process the set of observations as blocks (each
% observation could have 3 or less pieces of data)

if ~exist('gt','var')
    gt = [];
end

if ~weighted
    Lev_Fun = @(R) Fun( rot_input, R );
    [ R, err, errNorm, W ] = LM_Man_optim(Lev_Fun,R0,'space','SO(3)','debug',2, 'maxIters', 200);

else
    Lev_Fun = @(R) FunW( rot_input, R );
    [ R, err, errNorm, W ] = LM_Man_optim(Lev_Fun,R0,'space','SO(3)','weighted',true,'debug',2, 'maxIters', 200);   
end

figure, hold on
title(sprintf('Optimization with WEIGHTED=%d',weighted))
all_label( all_label==0 ) = [];
rgb = 'rgb';
subplot(211), hold on
for i=1:length(all_label)
    bar(i,err(i), rgb(all_label(i)))
end
subplot(212), hold on
for i=1:length(all_label)
    bar(i,W(i,i), rgb(all_label(i)))
end

end

function [residual, J] = Fun( rot_input, R )

n_cam = [rot_input.N];
line_LRF = [rot_input.l];
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
            MonteCarlo_simulate( )
        end
    end
end

% Build complete matrix from diagonal blocks
W = blkdiag(W{:});
% % Temporal debug:
% W = diag(diag(W));

    function MonteCarlo_simulate( )
        
        X0 = [ N(:); L(:) ];
        function e = F_e(X0)
            N = reshape(X0(1:3*3),3,3);
            L = reshape(X0(3*3+1:end),2,3);
            e = dot( N, R(1:3,1:2)*L, 1 )';
        end
        
        [~,J_man] = MonteCarlo.manDiffRot( N );
        A_N_man = J_man * A_N * J_man';
        % TODO: Input covariance A_N should be minimal covariance
        % NOTE: R for differential is not R_c_s but R_c_w=N -> Corrected
        fprintf('Montecarlo analysis of error covariance from N and L\n')
        [manSum,manSubtraction,manMean] = MonteCarlo.spaceOps( {'SO(3)','S1','S1','S1'}, {'R','R','R'} );
        tic, [mu,A] = MonteCarlo.simulate( @F_e, X0, blkdiag(A_N_man,A_L), 'manSum',manSum, 'manSubtraction',manSubtraction, 'manMean',manMean,...
            'AF', A_k, 'N', 1e4 ); toc
        keyboard
        
    end

end