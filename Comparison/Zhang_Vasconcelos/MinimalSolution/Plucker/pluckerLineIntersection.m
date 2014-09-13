% [i1, i2, d] = pluckerLineIntersection(L1, L2)
%
% Determines the intersection of two plucker lines, or its aproximation if 
% they do not intersect.
%
% USES:
%   [i1, i2, d] = pluckerLineIntersection(L1, L2)
%   [i, d]      = pluckerLineIntersection(L1, L2)
%   i           = pluckerLineIntersection(L1, L2)
%
% INPUT:
%   L1 - 6x1 vector with plucker coordinates of line 1
%   L2 - 6x1 vector with plucker coordinates of line 2
%
% OUTPUT:
%   i1 - Point from L1 with minimum distance to L2 
%   i2 - Point from L2 with minimum distance to L1 
%   i  - intersection of L1 and L2 (or its aproximation if L1 and L2 do not
%        intersect)
%   d  - distance between L1 and L2
%

function [i1, i2, d] = pluckerLineIntersection(L1, L2)

u1 = L1(1:3);
v1 = L1(4:6);
u2 = L2(1:3);
v2 = L2(4:6);

if u1.'*u1 ~= 0
    c1 = cross(v1,u1)/(u1.'*u1);
else
    c1 = [0;0;0];
end

if u2.'*u2 ~= 0
    c2 = cross(v2,u2)/(u2.'*u2);
else
    c2 = [0;0;0];
end

A = [u1 -u2];
b = c2-c1;
s = A\b;

i1 = c1 + s(1)*u1;
i2 = c2 + s(2)*u2;
d  = norm(i1-i2);

if nargout < 3
    i1 = (i1+i2)/2;
    i2 = d;
end