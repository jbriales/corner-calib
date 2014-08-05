function mu = manMeanS1( X )

mu = sum( X, 2 );
mu = mu / norm( mu );

end