function [residual, J] = f_Corner_2D( t, R, L, q, K )
% Err = f_Jac_Corner_2D( X, Rt )
% Implement here the computation of err vector (Nx1 for N elements in LSM
% problem)

N = size(L,2);
J = zeros(N,3);
residual = zeros(N,1);

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
    c = R(:,1:2) * c_3D + t;
    % Project the 3D point into the image plane
    c = K * c;
    % Calculate the Jacobian vector for each observation
    J(i,1) = fx*l(1)/(delta*c(3));
    J(i,2) = fy*l(2)/(delta*c(3));
    J(i,3) = ( c(3)*(cx*l(1)+cy*l(2)+l(3))-l'*c )/(delta*c(3)^2);
    % Estimates the residue
    residual(i) = (l'*c) / ( sqrt(l(1)^2+l(2)^2) * c(3) );    
end

end