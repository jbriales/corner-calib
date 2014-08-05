function out = Householder( x, a, ei )
% H  = Householder( x, ei )
% Ha = Householder( x, a, ei )
%   H is the Householder matrix associated to x
%   a is a column vector which multiplies the H matrix H*a (efficient product)

if ~exist('ei','var')
    ei = 3;
end
n = length( x );
I = eye( n );

v = x - sign( x(1) ) * norm(x) * I(:,ei);
if ~exist('a','var')
    H = I - 2 * v * v' / (v'*v);
    out = H;
else
    Ha = a - 2*v*(v'*a)/(v'*v);
    out = Ha;
end