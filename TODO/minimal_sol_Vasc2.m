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
    V{k} = null( kron(v{k}',[0 0 1]) );
    N{k} = [ null(n{k}') , n{k} ];
end

% Compute intersection of several nullspaces
xy = null( [N{1} -N{2}] );
% inter_12 = 
x = xy(1:size(N{1},2),:);
N12 = N{1}*x;
xy = null( [N12 N{3}] );
x = xy(1:size(N12,2),:);
N123 = N12*x;