function [rho, OmU, D] = computeTrihedronNormals2( Nbp )
% [rho, OmU, D] = computeTrihedronNormals2( Nbp )
% Solve the system of equations for trihedron directions in reduced bases

% Obtain null spaces (Om) in which v lies (orthogonality condition)
Nbp = mat2cell(Nbp,3,[1 1 1]); % Reshape to cell array
Om = cell(1,3);
for k=1:3
    Om{k} = null( Nbp{k}' );
end

% Compute bilinear forms matrices
B = cell(3,3);
for k=1:3
    ij = setdiff(1:3,k); i=ij(1); j=ij(2);
    B{i,j} = Om{i}' * Om{j};
    B{j,i} = Om{j}' * Om{i};
end

% Diagonalize bilinear forms
U = cell(1,3);
[U{1},~,~] = svd( B{1,2} );
[U{2},~,~] = svd( B{2,3} );
[U{3},~,~] = svd( B{3,1} );
% Make each base dextrorotatory
for k=1:3
    U{k} = U{k} * diag( [1 det(U{k})] );
end

D = cell(1,3);
d = zeros(1,3);
for k=1:3
    ij = setdiff(1:3,k); i=ij(1); j=ij(2);
    D{k} = U{i}'*B{i,j}*U{j};

    signD11 = sign(D{k}(1,1));
    d(k) = -D{k}(2,2) * signD11;
end

d_m = sqrt( prod(d) );
rho = snormalize( [ d_m d_m d_m ; d ] );

if nargout > 1
    % Return basis changes
    OmU = cell(1,3);
    for k=1:3
        OmU{k} = Om{k} * U{k};
    end
end