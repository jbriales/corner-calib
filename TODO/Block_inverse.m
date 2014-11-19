% Semi-diagonal Matrix inverse
Nblocks = 10;
blocks = cell(1,Nblocks);
blocks_inv = cell(1,Nblocks);
for k=1:Nblocks;
    blocks{k} = rand(3);
    blocks_inv{k} = inv(blocks{k});
end
main = blkdiag( blocks{:} );
main_inv = blkdiag( blocks_inv{:} );
mix = rand(3,3*Nblocks);
add = rand(3,3);

Complete = [ blkdiag( blocks{:} ), mix' ; mix, add ];

% Complete inverse
tic
sol1 = inv(Complete);
toc

% Block inverse
A = main;
A_inv = main_inv;
B = mix';
C = mix;
D = add;
tic
sol2 = [ inv(A-B*inv(D)*C) , -A_inv*B*inv(D-C*A_inv*B) ;
  -inv(D)*C*inv(A-B*inv(D)*C) , inv(D-C*A_inv*B) ];
toc

% Block inverse in low dimensions
A = main;
A_inv = main_inv;
B = mix';
C = mix;
D = add;
tic
perturbation = A_inv - A_inv * B * inv( -D + C * A_inv * B ) * C *A_inv;
sol3 = [ inv(A-B*(D\C)) , -A_inv*B*inv(D-C*A_inv*B) ;
  -inv(D)*C*perturbation , inv(D-C*A_inv*B) ];
toc

% In Matlab the computation time does not change (but maybe in other
% implementations it does)