% function [U,V,C,d]=GetLineInfo(L)
%
% GetLineInfo computes relevant information of line L (plucker
% coordinates). U is the direction versor, V is the unitary normal 
% of the plane cntaining L. C is the closest 3D point to the origin and 
% d is the distance

function [U,V,C,d]=GetLineInfo(L)

U=L(1:3);
V=L(4:6);
C=cross(V,U)/(transpose(U)*U);
u=norm(U);
v=norm(V);
d=u/v;
U=U*u^-1;
V=V*v^-1;
