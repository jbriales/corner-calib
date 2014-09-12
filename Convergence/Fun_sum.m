function [residual, J] = Fun_sum( corresp, x, K )

% TODO: vectorize the code... but i still don't know how 'kwak_input' is
% going to be

% Calculus of J(x) and E(x) for n matches
n = size(corresp,3); 
%K = Rig.Camera.K;
%T = transformation_expmap(x);   
T = [x; [0 0 0 1] ];
I = eye(3);
Z = zeros(3,3);
e1 = [1 0 0]';
e2 = [0 1 0]';
e3 = [0 0 1]';

B = kron(T',I);
C = [Z -skew(e1);
     Z -skew(e2);
     Z -skew(e3);
     I  Z];

D = B*C;

J = zeros(1,6);
residual = 0;

for i=1:n,
    
    % Homogeneous 3D points in the lidar frame
    p_s_c_h = makehomogeneous(corresp{2,1,i});
    p_s_l_h = makehomogeneous(corresp{2,2,i});
    p_s_r_h = makehomogeneous(corresp{2,3,i});
    
    
    % Transform the 3D points in the lidar frame to the image plane
    p_c_c = K * T(1:3,1:4) * p_s_c_h;
    p_c_l = K * T(1:3,1:4) * p_s_l_h;
    p_c_r = K * T(1:3,1:4) * p_s_r_h;
    
    % Estimate the distance to the line in the image plane (px)
    L_c = corresp{1,1,i};
    L_l = corresp{1,2,i};
    L_r = corresp{1,3,i};
    
    d_c = p_c_c' * L_c / p_c_c(3);
    d_l = p_c_l' * L_l / p_c_l(3);    
    d_r = p_c_r' * L_r / p_c_r(3);
    
    % Calculus of the residual
    residual = d_c^2 + d_l^2 + d_r^2 + residual;
    
    % Calculus of the jacobian
    kroen_c = kron(p_s_c_h', I);
    kroen_l = kron(p_s_l_h', I);
    kroen_r = kron(p_s_r_h', I);
    
    num_c   = K * T(1:3,1:4) * p_s_c_h ;    den_c   = num_c(3) ;
    num_l   = K * T(1:3,1:4) * p_s_l_h ;    den_l   = num_l(3) ;    
    num_r   = K * T(1:3,1:4) * p_s_r_h ;    den_r   = num_r(3) ;
    
    A_c     = L_c' * K * kroen_c / (den_c) - dot( L_c , num_c ) * ( K(3,:) * kroen_c ) / (den_c^2);
    A_l     = L_l' * K * kroen_l / (den_l) - dot( L_l , num_l ) * ( K(3,:) * kroen_l ) / (den_l^2);    
    A_r     = L_r' * K * kroen_r / (den_r) - dot( L_r , num_r ) * ( K(3,:) * kroen_r ) / (den_r^2);
    
    J       = J + d_c * A_c * D + d_l * A_l * D + d_r * A_r * D ;

end

J = 2 * J;

end
