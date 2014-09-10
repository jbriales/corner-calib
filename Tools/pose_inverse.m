function T_inv = pose_inverse(T)

R = T(1:3,1:3);  R_inv = inv(R);
T_inv = [ R_inv -R_inv*T(1:3,4); [ 0 0 0 1] ];

end