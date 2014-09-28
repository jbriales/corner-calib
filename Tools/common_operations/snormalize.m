function V = snormalize( V )
% Normalize columns of V to have norm 1 (lie on sphere)
N = size(V,1);
V = V ./ repmat( sqrt( sum( V.^2, 1 ) ), [N, 1, 1] ); % To use in 3D matrices too
end