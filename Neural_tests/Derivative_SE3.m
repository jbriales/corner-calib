% Script for testing SE(3) derivative
R90 = [ 0 -1; 1 0 ];

R = randn(3,3); [u,s,v] = svd(R); R = u*v'; clear u s v;
t = rand(3,1);
T_ = [R t];
T = [T_;0 0 0 1];

Os = zeros(3,3);
J = [Os, -skew(canvec(1))
     Os, -skew(canvec(2))
     Os, -skew(canvec(3))
     eye(3), Os];

kron(T',eye(3)) * J