function [A_l, A_p, A_ang, A_v, A_n] = propagateLineCov( pts, sigma )

% Uncertainty propagation
% For central point
N = size(pts,2);
A_p = sigma^2 * eye(2) / N;

% For direction
A_pts = kron( eye(N), sigma^2*eye(2) );

p = mean(pts,2);
L = pts - repmat(p,1,N);
M = cov(pts');
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
Hess_C_nn = -n'*(M'*M)*n + v'*(M'*M)*v;
Jac_C_nM  = (M*n)' * kron( -v',eye(2) ) + ...
    (-M*v)'* kron( +n',eye(2) );
J = - Hess_C_nn \ ( Jac_C_nM * J_M_p );
% Apply derivative to uncertainty propagation
A_ang = J * A_pts * J';
A_v   = n * A_ang * n';
A_n   = v * A_ang * v';

% For homogeneous line
% TODO: Take correlation into account
J_l_n = [ eye(2) ; -p' ];
J_l_p = [ zeros(2,2) ; -n' ];
A_l = J_l_n * v * A_ang * v' * J_l_n' + ...
      J_l_p * A_p * J_l_p';
end