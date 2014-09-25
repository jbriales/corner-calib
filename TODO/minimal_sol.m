clear

m = snormalize( randn(2,3) );
n = snormalize( randn(3,3) );

X = cell(1,3); N = cell(1,3);
for k=1:3
X{k} = kron( m(:,k), n(:,k) );
N{k} = null(X{k}');
end

xy = null( [N{1} N{2}] );
x = xy(1:size(N{1},2),:);
N12 = N{1}*x;
keyboard
xy = null( [N12 N{3}] );
x = xy(1:size(N12,2),:);
N123 = N12*x;
% (*) Do not admit decoupling from N123
% Divide common nullspace
N_1 = N123(1:3,:);
N_2 = N123(4:6,:);
Om = N_1' * N_2;