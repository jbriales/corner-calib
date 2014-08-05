function mu = manMeanRot( X )

R  = reshape( X, 3,3, [] );
[U,~,V] = svd( sum( R, 3 ) );
mu = U*V';
mu = mu(:);

end
