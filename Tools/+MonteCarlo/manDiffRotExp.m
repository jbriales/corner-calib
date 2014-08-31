function J_R_eps = manDiffRotExp( X )
% Jacobian of complete SO(3) space wrt tangent space so(3)

if isvector( X )
    R  = reshape( X, 3,3, [] );
elseif all( size(X) == [3 3] )
    R = X;
else
    error('Wrong dimensions')
end

J_R_eps  = [ -skew(R(:,1))
             -skew(R(:,2))
             -skew(R(:,3)) ];

end