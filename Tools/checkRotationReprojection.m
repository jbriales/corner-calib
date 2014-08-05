function checkRotationReprojection( c, R_c_w, K, format )
% checkRotationReprojection( c, R_c_w )

if ~exist('format','var')
    format = '--';
end

L = mat2cell( skew(makehomogeneous(c)) * R_c_w, 3,[1 1 1]);
[c,L] = uncalibrateData(c,L,K);

% imshow( img.I ); hold on
plot(c(1),c(2),'om');

rgb = 'rgb';
for k=1:3
    plotHomLineWin( L{k}, [rgb(k),format] )
end