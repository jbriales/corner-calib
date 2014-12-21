function L = lnormalize( L )
% Normalize columns of L for normal vector to have norm 1 (lie on circle)
% This normalization only makes sense for homogeneous lines in P2
if size(L,1) ~= 3
    error('Columns must be 3-vectors');
end
L = L ./ repmat( sqrt( sum( L(1:2).^2, 1 ) ), [3, 1, 1] ); % To use in 3D matrices too
end