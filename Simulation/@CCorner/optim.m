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

% R_aux = Rig.R_c_s + randn(3,3)*0.01;
% [U,S,V] = svd(R_aux);
% x0 = [U*V' [0.15 0 0]'];

if ~weighted
    Lev_Fun = @(x) Fun( corresp, x, K );
    [ x, err, errNorm, W ] = LM_Man_optim(Lev_Fun,x0,'space','SE(3)','debug',2, 'maxIters', 200);
else
    Lev_Fun = @(x) FunW( corresp, x, K, theta );
    [ x, err, errNorm, W ] = LM_Man_optim(Lev_Fun,x0,'space','SE(3)','debug',2, 'maxIters', 200);   
end

end

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
    residual = d_l^2 + d_c^2 + d_r^2 + residual;
    
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
    
    J       = J + A_l * D + A_c * D + A_r * D ;

end

end

function [residual, J] = FunW( corresp, x, K, theta )

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

d   = zeros(3,n);
d_d = zeros(3,n,6);
d_e = zeros(3,n,6);
H   = zeros(3,n);
d_H = zeros(3,n);

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
    
    % Calculus of the weight
    d(1,i)   = d_l^2;
    d(2,i)   = d_c^2;
    d(3,i)   = d_r^2;
    
    % Calculus of the huber penalty function
    [H_aux, d_H_aux] = huber_penalty(d(:,i), d_max);
    H(:,i) = H_aux;
    d_H(:,i) = d_H_aux;
       
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
    
    for j = 1:6    
        d_d(1,i,j) = A_l(j);    d_e(1,i,j) = d(1) * d_d(1,i,j);
        d_d(2,i,j) = A_c(j);    d_e(2,i,j) = d(1) * d_d(2,i,j);
        d_d(3,i,j) = A_r(j);    d_e(3,i,j) = d(1) * d_d(3,i,j);
    end
    
end

% Calculus of the jacobian and the residual
e = sum(d.^2,2)/n;
w = 1./e;

d_e_s = sum(d_e,2); reshape(d_e_s,3,6);

for i = 1:n
    J_ = zeros(1,6); residual_ = 0;
    for j = 1:3
        
        J_ =  J_ -w(j)^2 * d_e_s(j,:) * H(j,i) / n + w(j) * d_H(j,i) * reshape(d_d(j,i,:),1,6);
        residual_ = residual_ + w(j) * H(j,i);        
        
    end
    J(i,:) = J_;
    residual(i,1) = residual_;
    
end

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

d   = zeros(3,n);
d_d = zeros(3,n,6);
d_e = zeros(3,n,6);
H   = zeros(3,n);
d_H = zeros(3,n);

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
    
    % Calculus of the weight
    d(1,i)   = d_l^2;
    d(2,i)   = d_c^2;
    d(3,i)   = d_r^2;
    
    % Calculus of the huber penalty function
    [H_aux, d_H_aux] = huber_penalty(d(:,i), d_max);
    H(:,i) = H_aux;
    d_H(:,i) = d_H_aux;
       
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
    
    for j = 1:6    
        d_d(1,i,j) = A_l(j);    d_e(1,i,j) = d(1) * d_d(1,i,j);
        d_d(2,i,j) = A_c(j);    d_e(2,i,j) = d(1) * d_d(2,i,j);
        d_d(3,i,j) = A_r(j);    d_e(3,i,j) = d(1) * d_d(3,i,j);
    end
    
end

% Calculus of the jacobian and the residual
e = sum(d.^2,2)/n;
w = 1./e;

d_e_s = sum(d_e,2); reshape(d_e_s,3,6);

for i = 1:n
    for j = 1:3
        
        J = J -w(j)^2 * d_e_s(j,:) * H(j,i) / n + w(j) * d_H(j,i) * reshape(d_d(j,i,:),1,6);
        residual = residual + w(j) * H(j,i);        
        
    end
end

end

function [H, d_H]  =  huber_penalty(d, d_max)

for i = 1:3
    
    if d(i) >= d_max
        H(i,1)   = d(i) ^ 2 ;
        d_H(i,1) = 2 * d(i) ;
    else
        H(i,1)   = d_max * ( 2 * abs(d(i)) - d_max ) ;
        d_H(i,1) = 2 * d_max * sign(d(i)) ;  
    end
    
end

end