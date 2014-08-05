function [residual, J] = f_Corner_3D( t, R, L, q )
% Err = f_Jac_Corner_3D( X, Rt )
% Implement here the computation of err vector (Nx1 for N elements in LSM
% problem)

N = size(L,2);
J = zeros(N,3);
residual = zeros(N,1);

R = R(1:3,1:2);

for i=1:N
    l = L(:,i);
    %l = K' \ l;
    c_3D = q(:,i);
    alpha = norm(l(1:2));
    
    % Transform c_3D into Camera frame (equivalent to homogeneous image point)
    c_w = R(:,1:2) * c_3D + t;
    % Calculate the Jacobian vector for each observation
    J(i,1) = l(1)/alpha;
    J(i,2) = l(2)/alpha;
    J(i,3) = l(3)/alpha;
    % Estimates the residue
    residual(i) = (l'*c_w) / alpha;        
end

end