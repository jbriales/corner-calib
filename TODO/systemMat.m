% List of monomials
% x3z, xy2z, xz3
% x3, x2y, xy2, y3, x2z, xyz, y2z, xz2, yz2, z3
% x2, xy, y2, xz, yz, z2
% x, y, z
% 1
% M = ...
%     [ zeros(1,3) , zeros(1,10) , D(1,1), 0, D(1,2), 0, 0, D(1,3), zeros(1,3) -1 ;
%       zeros(1,3) , zeros(1,10) , D(2,1), 0, D(2,2), 0, 0, D(2,3), ;

B = 0.5*(Bd+Bd');
vB = [B(1,2), B(2,3), B(3,1)];
D = [ diag(Da), diag(Db), diag(B) ]';

M = [ zeros(1,7), D(1,1) ;
      zeros(1,7), D(1,2) ;
      zeros(1,7), D(1,3) ;
      %
      zeros(1,3), D(3,1), D(1,1), zeros(1,3) ;
      zeros(1,3), 2*B(1,2), D(1,1), 0, 0, 0 ;
      zeros(1,3), D(3,2), D(1,2), 0, 0, 0 ;
      zeros(1,5), D(1,2), zeros(1,2) ;
      zeros(1,3), 2*B(1,3), 0 0 , D(1,1), 0 ;
      zeros(1,3), 2*B(2,3), zeros(1,4) ;
      zeros(1,3), 0 0 0 , D(1,2), 0 ;
      zeros(1,3), D(3,3), D(1,3), 0 0 0 ;
      zeros(1,3), 0 0 , D(1,3), 0 0 ;
      zeros(1,3), 0 0 0 D(1,3) 0 ;
      %
      D(:,1)' zeros(1,5) ;
      0 0 2*B(1,2) 0 0 0 0 0 ;
      D(:,2)' zeros(1,5) ;
      0 0 2*B(1,3) 0 0 0 0 -1 ;
      0 0 2*B(2,3) 0 0 0 0 0 ;
      D(:,3)' zeros(1,5) ;
      %
      zeros(1,4) -1 0 0 0 ;
      zeros(1,4) 0 -1 0 0 ;
      zeros(1,4) 0 0 -1 0 ;
      %
      -1 -1 zeros(1,6) ]';
  
  % GT value: vector of monomials
  x=rho(1);y=rho(2);z=rho(3);
  vec_mons = [ x^3*z
               x*y^2*z
               x*z^3
               x^3
               x^2*y
               x*y^2
               y^3
               x^2*z
               x*y*z
               y^2*z
               x*z^2
               y*z^2
               z^3
               x^2
               x*y
               y^2
               x*z
               y*z
               z^2
               x
               y
               z
               1 ];
      
% Permutation matrix
new = 1:23;
old = [1 2 4 5 6 7 10 12 16 , 3 8 9 11 14 15 17 20 , 13 18 19 21 22 23 ];
T = zeros(23,23);
for i=1:23
    T(old(i),new(i)) = 1;
end

% Reordered matrix
M_ = M * T;

% Solve with built-in function
syms x y z
% rho = [x; y; z];
% x = rho(1), y = rho(2), z = rho(3)
rho2 = [x^2; y^2; z^2];
rhoxy = [x*y; y*z; z*x];
eq1 = D(1,:) * rho2 == 1;
eq2 = D(2,:) * rho2 == 1;
eq3 = D(3,:) * rho2 + 2*vB * rhoxy == 0;
% S = solve( eq1, eq2, eq3, x,y,z )
[x,y,z] = solve( eq1, eq2, eq3, x,y,z );
sols = [x y z].';
sols2 = sols.^2;