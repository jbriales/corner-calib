function [residual, J, W] = f_Corner_2D_W( t, R, L, q, K, A_lh, A_q)
% Err = f_Jac_Corner_2D( X, Rt )
% Implement here the computation of err vector (Nx1 for N elements in LSM
% problem)

N = size(L,2);
J = zeros(N,3);
W = zeros(N,N);
residual = zeros(N,1);

R_full = R;
R = R(1:3,1:2);
% t = Rt(1:3,4);

fx = K(1,1);
fy = K(2,2);
cx = K(1,3);
cy = K(2,3);

for i=1:N
    l = L(:,i);
    l = K' \ l;
    c_3D = q(:,i);
    delta = norm(l(1:2));
    
    % Transform c_3D into Camera frame (equivalent to homogeneous image point)
    c_w = R(:,1:2) * c_3D + t;
    % Project the 3D point into the image plane
    c = K * c_w;
    % Calculate the Jacobian vector for each observation
    J(i,1) = fx*l(1)/(delta*c(3));
    J(i,2) = fy*l(2)/(delta*c(3));
    J(i,3) = ( c(3)*(cx*l(1)+cy*l(2)+l(3))-l'*c )/(delta*c(3)^2);
    % Estimates the residue
    residual(i) = (l'*c) / ( sqrt(l(1)^2+l(2)^2) * c(3) );    
    
    % Estimates the weights based on the uncertainties    
%     J_lh_l = [0 -1; 1 0; -p(2,i) p(1,i)];
%     J_lh_p = [0 0 ; 0 0; lh(2,i) lh(1,i)];
%     J_l = [0; 1];
%     A_l_ = J_l * A_l{i} * J_l';
%    
%     A_lh = J_lh_l* A_l_ *J_lh_l' + J_lh_p* A_p{i} *J_lh_p';
%     
%     
%     A_q = eye(2);   %TODO

    
    J_qc_qs = [fx/c_w(3) 0 -fx*c_w(1)/c_w(3)^2;
               0 fy/c_w(3) -fy*c_w(2)/c_w(3)^2]*R_full;
    A_qc = J_qc_qs * A_q{i} * J_qc_qs';
    
    J_d_qc = [l(1) l(2)];
    J_d_lh  = [c(1) c(2) 1];
    A_d = J_d_qc * A_qc * J_d_qc' + J_d_lh * A_lh{i} * J_d_lh';
    
    W(i,i) = 1/A_d;    
    
 end

end