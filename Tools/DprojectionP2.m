function J = DprojectionP2( p )
% J = DprojectionP2( p )
% Derivative of the P2-R2 projection u = x/z, v = y/z, very common when
% working with images
% Input: 3x1 homogeneous vector
% Output: 2x3 jacobian
if numel(p)==2
    p = makehomogeneous(p);
end

J = 1/p(3)^2 * [ p(3) 0 -p(1) ;
                 0 p(3) -p(2) ];