function inc = manSubtractionRn( X, X0 )

N   = size(X,2);
X0  = repmat(X0,1,N);
inc = X-X0;

end