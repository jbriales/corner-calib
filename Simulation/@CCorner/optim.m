function [ x, err, errNorm, W ] = optim( obj, corresp, x0, weighted, Rig )

% :: Inputs
%
%   corresp  - cell array ( 2, 3, N ) with the input data (PRE-PROCESS)
%   x0       - initial estimation of the SE(3) transformation (3x4 matrix)
%   weighted - parameter to choose the optimization method
%
% :: Outputs
%   x        - optimal estimation of the SE(3) transformation (6-vector)

K = Rig.Camera.K;
theta = Rig.Lidar.FOVd / Rig.Lidar.N ;

% R0 = [ 0 -1  0
%        0  0 -1
%        1  0  0 ];
% x0 = [R0 zeros(3,1)];
x0 = [Rig.R_c_s Rig.t_c_s];
% x0 = [Rig.R_c_s [0.16 0.001 -0.001]'];
 
R_aux = Rig.R_c_s + randn(3,3)*0.01;
[U,S,V] = svd(R_aux);
x0 = [U*V' [0.15 0 0]'];

if ~weighted
    Lev_Fun = @(x) Fun( corresp, x, K );
    [ x, err, errNorm, W ] = LM_Man_optim(Lev_Fun,x0,'space','SE(3)','debug',0, 'maxIters', 200);
else
    Lev_Fun = @(x) FunW( corresp, x, K, theta );
    [ x, err, errNorm, W ] = LM_Man_optim(Lev_Fun,x0,'space','SE(3)','weighted',true,'debug',0, 'maxIters', 200);   
end

end

function [residual, J] = Fun( corresp, x, K )

% TODO: vectorize the code... but i still don't know how 'kwak_input' is
% going to be

% Calculus of J(x) and E(x) for n matches
n = size(corresp,3); 
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

J        = zeros(3*n,6);
residual = zeros(3*n,1);

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
    residual(3*i-2,1) = d_c;
    residual(3*i-1,1) = d_l;
    residual(3*i,1)   = d_r;
    
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
          
    J(3*i-2,:) = 2 * d_c * A_c * D;
    J(3*i-1,:) = 2 * d_l * A_l * D;
    J(3*i,:)   = 2 * d_r * A_r * D;
    

end

end

function [residual, J, W] = FunW( corresp, x, K, theta )

% TODO: vectorize the code... but i still don't know how 'kwak_input' is
% going to be
% TODO: their weight is not like ours, so the optimization is the same but
% the residuals change

% Calculus of J(x), E(x) and W(x) for n matches
n = size(corresp,3);
% K = Rig.Camera.K;
% T = transformation_expmap(x);
T = [x; [0 0 0 1] ];
I = eye(3);
Z = zeros(3,3);
e1 = [1 0 0]';
e2 = [0 1 0]';
e3 = [0 0 1]';
%theta   = Rig.Lidar.FOVd / Rig.Lidar.N ; 
d_max = K(1,1) * tand(theta/2);

B = kron(T',I);
C = [Z -skew(e1);
     Z -skew(e2);
     Z -skew(e3);
     I  Z];

D = B*C;

J        = zeros(3*n,6);
residual = zeros(3*n,1);
W        = zeros(3*n,3*n);

omega = zeros(3,1);
H_n   = zeros(3,1);

alpha = zeros(3*n,6);
gamma = zeros(3*n,6);

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
   
    % Calculus of the weights
    omega(1,1)   = omega(1,1) + d_c^2;
    omega(2,1)   = omega(2,1) + d_l^2;
    omega(3,1)   = omega(3,1) + d_r^2;
    
    % Calculus of the huber penalty function
    [H, d_H] = huber_penalty([d_c d_l d_r]' , d_max);
    H_n      = H_n + H;
    
    % Calculus of the residual (TODO)
    residual(3*i-2,1) = sign(H(1))*sqrt(abs(H(1)));
    residual(3*i-1,1) = sign(H(2))*sqrt(abs(H(2)));
    residual(3*i,1)   = sign(H(3))*sqrt(abs(H(3)));
%     residual(3*i-2,1) = d_c;
%     residual(3*i-1,1) = d_l;
%     residual(3*i,1)   = d_r;
           
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
    
    alpha(3*i-2,:) = d_H(1) * A_c * D;
    alpha(3*i-1,:) = d_H(2) * A_l * D;
    alpha(3*i,:)   = d_H(3) * A_r * D;
    
    gamma(3*i-2,:) = d_c * A_c * D ;
    gamma(3*i-1,:) = d_l * A_l * D ;
    gamma(3*i,:)   = d_r * A_r * D ;  
  
end

% Calculus of the weight function
omega = n./omega;

% Calculus of the jacobian and the weight matrices
for i = 1:n
    % Full jacobian
    J(3*i-2,:) = omega(1) * alpha(3*i-2,:) - 2 * omega(1)^2 * H_n(1) * gamma(3*i-2,:) / n;
    J(3*i-1,:) = omega(2) * alpha(3*i-1,:) - 2 * omega(2)^2 * H_n(2) * gamma(3*i-1,:) / n;
    J(3*i,:)   = omega(3) * alpha(3*i,:)   - 2 * omega(3)^2 * H_n(3) * gamma(3*i,:) / n;
      
%     % Simplified jacobian
%     J(3*i-2,:) = 2 * omega(1) * gamma(3*i-2,:);
%     J(3*i-1,:) = 2 * omega(2) * gamma(3*i-1,:);
%     J(3*i,:)   = 2 * omega(3) * gamma(3*i,:);
    
    W(3*i-2,3*i-2) = omega(1);
    W(3*i-1,3*i-1) = omega(2);
    W(3*i,3*i)     = omega(3);

end



end




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
% residual = sign(residual) * sqrt(abs(residual));

end

function [residual, J] = FunW_sum( corresp, x, K, theta )

% TODO: vectorize the code... but i still don't know how 'kwak_input' is
% going to be
% TODO: their weight is not like ours, so the optimization is the same but
% the residuals change

% Calculus of J(x), E(x) and W(x) for n matches
n = size(corresp,3);
% K = Rig.Camera.K;
% T = transformation_expmap(x);
T = [x; [0 0 0 1] ];
I = eye(3);
Z = zeros(3,3);
e1 = [1 0 0]';
e2 = [0 1 0]';
e3 = [0 0 1]';
%theta   = Rig.Lidar.FOVd / Rig.Lidar.N ; 
d_max = K(1,1) * tand(theta/2);

B = kron(T',I);
C = [Z -skew(e1);
     Z -skew(e2);
     Z -skew(e3);
     I  Z];

D = B*C;

J = zeros(1,6);
residual = 0;

omega = zeros(3,1);
H_n   = zeros(3,1);

alpha = zeros(3,6);
gamma = zeros(3,6);

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
    
    % Calculus of the weights
    omega(1,1)   = omega(1,1) + d_c^2;
    omega(2,1)   = omega(2,1) + d_l^2;
    omega(3,1)   = omega(3,1) + d_r^2;
    
    % Calculus of the huber penalty function
    [H, d_H] = huber_penalty([d_c d_l d_r]' , d_max);
    H_n      = H_n + H;
           
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
    
    alpha(1,:) = alpha(1,:) + d_H(1) * A_c * D;
    alpha(2,:) = alpha(2,:) + d_H(2) * A_l * D;
    alpha(3,:) = alpha(3,:) + d_H(3) * A_r * D;
    
    gamma(1,:) = gamma(1,:) + d_c * A_c * D ;
    gamma(2,:) = gamma(2,:) + d_l * A_l * D ;
    gamma(3,:) = gamma(3,:) + d_r * A_r * D ;  
  
end

% Calculus of the weight function
omega = n./omega;

% Calculus of the jacobian and the residual
residual = omega' * H_n;
for k = 1:3
    %J = J + omega(k) * alpha(k,:) - 2 * omega(k)^2 * H_n(k) * gamma(k,:) / n;   
	J = J + omega(k) * alpha(k,:);   
end

% residual = sign(residual) * sqrt(abs(residual));

end

function [H, d_H]  =  huber_penalty(d, d_max)

for i = 1:3
    
    if abs(d(i)) >= d_max
        H(i,1)   = d(i) ^ 2 ;
        d_H(i,1) = 2 * d(i) ;
    else
        H(i,1)   = d_max * ( 2 * abs(d(i)) - d_max ) ;
        d_H(i,1) = 2 * d_max * sign(d(i)) ;  
    end
    
end

end
