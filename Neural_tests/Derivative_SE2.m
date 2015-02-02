% Script for testing SE(2) derivative
R90 = [ 0 -1; 1 0 ];

n = snormalize(rand(2,1));
n_ = R90 * n;
R = [n , n_];
t = rand(2,1);
T_ = [R t];
T = [T_;0 0 1];

Je = [0 0 0 ; 0 0 1 ; 0 0 -1 ; 0 0 0 ; 1 0 0 ; 0 1 0 ];
Jl = ( Je * diag([1 1 0.5]) )';

J = kron(T',eye(2)) * Je
iJ = Jl * kron(inv(T'),eye(2))

R90R = R90 * R;
iJ = [ kron((-R'*t)',eye(2)), eye(2)
       0.5* R90R(:)', 0 0 ]
   
iJ = [ kron((-R'*t),eye(2)), 0.5* R90R(:);
       eye(2), zeros(2,1) ]'