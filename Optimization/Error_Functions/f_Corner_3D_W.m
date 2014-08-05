function [residual, J, W] = f_Corner_3D_W( t, R, L, q, A_lh, A_q)
% Err = f_Jac_Corner_3D( X, Rt )
% Implement here the computation of err vector (Nx1 for N elements in LSM
% problem)

N = size(L,2);
J = zeros(N,3);
W = zeros(N,N);
residual = zeros(N,1);

R = R(1:3,1:2);

for i=1:N
    l = L(:,i);
    c_3D = q(:,i);
    alpha = norm(l); alpha2 = alpha^2; alpha3 = alpha^3;
    
    % Transform c_3D into Camera frame (equivalent to homogeneous image point)
    c_w = R(:,1:2) * c_3D + t;
    % Calculate the Jacobian vector for each observation
    J(i,1) = l(1)/alpha;
    J(i,2) = l(2)/alpha;
    J(i,3) = l(3)/alpha;
    % Estimates the residue
    r = (l'*c_w);
    residual(i) = r / alpha;    
    
    % Estimates the weights based on the uncertainties  
    J_d_q  = J(i,:);
    J_d_lh  = [c_w(1)*alpha2-r*l(1) c_w(2)*alpha2-r*l(2) c_w(3)*alpha2-r*l(3) ] / alpha3;
    A_d = J_d_q * A_q{i} * J_d_q' + J_d_lh * A_lh{i} * J_d_lh';
    W(i,i) = 1/A_d;    
    
 end

end