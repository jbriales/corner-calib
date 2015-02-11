function Hn = Reflection( n )
% Hn = Reflection( n )
% Return the symmetric matrix I-2nn'/n'n (see Householder matrix)
% which reflects any vector x wrt the plane whose normal vector is n

if any( size(n) ~= [3 1] )
    error('Bad size for argument n, should be 3x1');
end
Hn = eye(3) - 2 * (n * n') / (n'*n);
end