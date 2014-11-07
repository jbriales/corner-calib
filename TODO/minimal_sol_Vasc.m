% Minimal solution (alternative to Vasconcelos)
clear
clc

% Data is generated through simulation with GT
load(fullfile(pwd,'TODO','data.mat'), 'n','v','R');
% Check solution
% dot(n,R(:,1:2)*v,1)

n = num2cell(n(:,1:3),1);
v = num2cell(v(:,1:3),1);
% n = num2cell( snormalize(rand(3)),1 );
% v = num2cell( snormalize(rand(2,3)),1 );

% for k=1:3
%     V{k} = [ null(n{k}') n{k} ];
% end

% Linear equation for every diagonalized representation
for k=1:3
    eq_data{k} = kron( v{k}', [0 0 1] );
end

% Linear equation for equality of r's
% for k=1:3
%     ij = num2cell(setdiff(1:3,k));
%     [i,j] = deal( ij{:} );
%     eq_repr{k} = [ V{i} -V{j} ];
% end

% Get bi-orthogonal nullspace
for k=1:3
    N{k} = null( kron( v{k}', [0 0 1] ) );
end

% Compute intersection of several nullspaces
xy = null( [N{1} N{2}] );
x = xy(1:size(N{1},2),:);
N12 = N{1}*x;
xy = null( [N12 N{3}] );
x = xy(1:size(N12,2),:);
N123 = N12*x;

% Convert N to bi-orthogonal
NN = N;
for k=1:3
    b = [N{k}(end,2) N{k}(end,end)]';
    A = [N{k}(1:3,2) N{k}(1:3,end)];
    [u,s,~] = svd( A'*A );
    
    ub = u*b;
    psi_y = 1/ub(2);
    psi_x = sqrt(1-s(end)*psi_y^2);
    psi1 = [+psi_x psi_y]';
    psi2 = [-psi_x psi_y]';
    % delta = u*psi;
    eps1 = A * u * psi1;
    eps2 = A * u * psi2;
    
    % Transform nullspace to new bi-orthogonal basis
    NN{k}(:,2)   = [eps1' 0 0 1]';
    NN{k}(:,end) = [eps2' 0 0 1]';
end