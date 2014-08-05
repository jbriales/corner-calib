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
n_cam = [rot_input.N];
line_LRF = [rot_input.l];
residual = dot( n_cam, R(1:3,1:2) * line_LRF, 1 )';
J = cross( R(1:3,1:2) * line_LRF, n_cam, 1 )';

% Weights
N_obs = length(rot_input);
W = cell(1, N_obs);

ort = [ 0 -1
        1  0 ];
    
R = R(1:3,1:2); % Only 2 columns are used

for i=1:N_obs
    N = rot_input(i).N;
    A_N = rot_input(i).A_N;
    L = rot_input(i).l;
    A_L = rot_input(i).A_l;
    if size( N, 2 ) == 3
%         i = 3*(k-1);
%         span = (1:3)+i;
%         iN = 9*(k-1);
%         spanN = (1:9)+iN;
        n1 = N(:,1);
        n2 = N(:,2);
        n3 = N(:,3);
        l1 = L(:,1);
        l2 = L(:,2);
        l3 = L(:,3);
        
        JN  = blkdiag( l1'*R', l2'*R', l3'*R' );
        Ja1 = n1' * R * ort * l1;
        Ja2 = n2' * R * ort * l2;
        Ja3 = n3' * R * ort * l3;
        
        A_Nk = JN * A_N * JN';
        A_Lk = diag( [Ja1 Ja2 Ja3].^2 ) * A_L;
        A_k  = A_Nk + A_Lk;
    elseif size( N, 2 ) == 2
        n1 = N(:,1);
        n2 = N(:,2);
        l1 = L(:,1);
        l2 = L(:,2);
        
        JN  = blkdiag( l1'*R', l2'*R' );
        Ja1 = n1' * R * ort * l1;
        Ja2 = n2' * R * ort * l2;
        
        A_Nk = JN * A_N * JN';
        A_Lk = diag( [Ja1 Ja2].^2 ) * A_L;
        A_k  = A_Nk + A_Lk;
    else
        n1 = N(:,1);
        l1 = L(:,1);
        
        JN  = blkdiag(l1'*R');
        Ja1 = n1' * R * ort * l1;
        
        A_Nk = JN * A_N * JN';
        A_Lk = diag( Ja1^2 ) * A_L;
        A_k  = A_Nk + A_Lk;
    end
    W{i} = pinv(A_k);
end

% Build complete matrix from diagonal blocks
W = blkdiag(W{:});
% % Temporal debug:
% W = diag(diag(W));

end