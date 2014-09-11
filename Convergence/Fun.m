
function [residual, J] = Fun( corresp, x, K )

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

J = zeros(n,6);
residual = zeros(n,1);

for i=1:n,
    
    % Transform the 3D points in the lidar frame to the image plane
    p_c_c = K * T(1:3,1:4) * makehomogeneous(corresp{2,1,i});    p_c_c_h = makehomogeneous( p_c_c );
    p_c_l = K * T(1:3,1:4) * makehomogeneous(corresp{2,2,i});    p_c_l_h = makehomogeneous( p_c_l );
    p_c_r = K * T(1:3,1:4) * makehomogeneous(corresp{2,3,i});    p_c_r_h = makehomogeneous( p_c_r );
    
    % Estimate the distance to the line in the image plane (px)   TODO
    d_c = p_c_c' * corresp{1,1,i} / p_c_c(3);
    d_l = p_c_l' * corresp{1,2,i} / p_c_l(3);    
    d_r = p_c_r' * corresp{1,3,i} / p_c_r(3);
    
    L_c = corresp{1,1,i};
    L_l = corresp{1,2,i};
    L_r = corresp{1,3,i};
    
    % Calculus of the residual
    residual(i,1) = d_l^2 + d_c^2 + d_r^2 ;
    
    % Calculus of the jacobian
    kroen_l = kron(p_c_l_h', I);
    kroen_c = kron(p_c_c_h', I);
    kroen_r = kron(p_c_r_h', I);
    
    num_l   = K * T(1:3,1:4) * p_c_l_h ;    den_l   = num_l(3) ;
    num_c   = K * T(1:3,1:4) * p_c_c_h ;    den_c   = num_c(3) ;
    num_r   = K * T(1:3,1:4) * p_c_r_h ;    den_r   = num_r(3) ;
    
    A_l     = L_l' * K * kroen_l / (den_l) - dot( L_l , num_l ) * ( K(3,:) * kroen_l ) / (den_l^2);
    A_c     = L_c' * K * kroen_c / (den_c) - dot( L_c , num_c ) * ( K(3,:) * kroen_c ) / (den_c^2);
    A_r     = L_r' * K * kroen_r / (den_r) - dot( L_r , num_r ) * ( K(3,:) * kroen_r ) / (den_r^2);
    
    J(i,:)  = A_l * D + A_c * D + A_r * D ;

end

end
