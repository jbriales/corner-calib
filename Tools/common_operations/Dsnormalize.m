function J = Dsnormalize( v )
% J = Dsnormalize( v )
% Compute derivative of vector normalization operation v/norm(v)
% IMPORTANT: The derivative should be always projected to manifold tangent space

N = size(v,1); % Vector should be column
J = 1/norm(v) * ( eye(N) - v*v'/(v'*v) );