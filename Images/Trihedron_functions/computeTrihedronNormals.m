function obj_Rtri = computeTrihedronNormals( xi, K, Nbp )
% R_tri = computeTrihedronNormals( c, K, v_im )
% 

Dproj = DprojectionP2( xi.c.X );
Ort   = [ 0 -1 ; 1 0 ];
v_im  = [ xi.v(1).X , xi.v(2).X , xi.v(3).X ];
% Obtain null spaces (Om) in which n lies (parallel condition)
for i=1:3
    Om{i} = null( (Ort * v_im(:,i))' * Dproj * K );
    Om{i} = null( (Ort * v_im(:,i))' * Dproj * K );
end

% Compute bilinear forms matrices
Q = cell(3,3);
q = cell(1,3);
for k=1:3
    temp = num2cell( setdiff(1:3,k) );
    [i,j] = deal( temp{:} );
    
    q{k}   = Om{i}' * Om{j};
    Q{k,k} = zeros(2,2);
    Q{i,j} = Om{i}'*Om{j};
    Q{j,i} = Om{j}'*Om{i};
end
% Diagonalize bilinear forms
U = cell(1,3);
[U{1},~,~] = svd( Q{1,2} );
[U{2},~,~] = svd( Q{2,3} );
[U{3},~,~] = svd( Q{3,1} );
% Make each base dextrorotatory
for k=1:3
    U{k} = U{k} * diag( [1 det(U{k})] );
end
%             [U_{1},~,~] = svd( Q{1,3} );
%             [U_{2},~,~] = svd( Q{2,1} );
%             [U_{3},~,~] = svd( Q{3,2} );
%             for i=1:3 % Check equal matrices
%                 if ~isequal(U{i},U_{i})
%                     warning('Different U decompositions');
%                 end
%             end
D = cell(1,3);
d = zeros(1,3);
for k=1:3
    temp = num2cell( setdiff(1:3,k) );
    [i,j] = deal( temp{:} );
    D{k} = U{i}'*Q{i,j}*U{j};
    % Is possible to change sign without more concern?
    D{k} = D{k} * sign(D{k}(1,1));
    d(k) = -D{k}(2,2);
end
%             keyboard
d_m = sqrt( prod(d) );
xp = +d_m ./ d;
xm = -d_m ./ d;
xpm = {xp, xm};
for i_sign = 1:2 % Compute with both signs and choose det = 1
    x = xpm{i_sign};
    
    y = ones(1,3);
    rho = [x ; y];
    rho = mat2cell( rho, 2,[1 1 1] );
    lambda = cell(1,3);
    for k=1:3
        lambda{k} = snormalize( U{k} * rho{k} );
    end
    % Compute sign of dot product, correct, and normals
    s = zeros(1,3);
    n = cell(1,3);
    for i=1:3
        s(i) = sign( v_im(:,i)' * Dproj * K * Om{i} * lambda{i} );
        lambda{i} = s(i) * lambda{i};
        n{i} = Om{i} * lambda{i};
    end
    R_tri = cell2mat( n );
    
    if abs( det(R_tri) - 1 ) < 0.5
        break
    end
end

if 0 % Temporal debug (check with other sign)
    x_ = xm;
    rho_ = [x_ ; y];
    rho_ = mat2cell( rho_, 2,[1 1 1] );
    lambda_ = cell(1,3);
    for k=1:3
        lambda_{k} = snormalize( U{k} * rho_{k} );
    end
    % Compute sign of dot product, correct, and normals
    s_ = zeros(1,3);
    n_ = cell(1,3);
    for i=1:3
        s_(i) = sign( v_im(:,i)' * Dproj * K * Om{i} * lambda_{i} );
        lambda_{i} = s_(i) * lambda_{i};
        n_{i} = Om{i} * lambda_{i};
    end
    R_tri_ = cell2mat( n_ );
end

% Compute covariance: Temporarily use old function Phi
obj_Rtri = Manifold.SO3( R_tri );

if exist('Nbp','var')
    % Temporarily use the other function
    N = Nbp.arr;
    J_Phi = @(R) cross( R, N, 1 )';
    % Compute covariance of trihedron normals
    J_Phi_eps = J_Phi( R_tri );
    J_Phi_Nbp = mat2cell( R_tri, 3, [1 1 1] );
    J_Phi_Nbp = blkdiag( J_Phi_Nbp{:} )';
    J_eps_Nbp = - J_Phi_eps \ J_Phi_Nbp;
    obj_Rtri.setMinimalCov( J_eps_Nbp * Nbp.A_X * J_eps_Nbp' );
end

end