function plotErrWeights( Err, W )

figure, hold on
bar(Err)

Err_max = max(Err);
iW = 1./diag(W);
diagWeights = diag(iW)/ max(diag(iW)) * 10*Err_max .* sign( Err );
plot( diagWeights, 'r*' )