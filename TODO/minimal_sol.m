clear

% Data is generated through simulation with GT
load(fullfile(pwd,'TODO','data.mat'), 'n','v','R');
% Check solution
% dot(n,R(:,1:2)*v,1)
r1 = R(:,1);
r2 = R(:,2);
r  = [r1' r2']';

n = num2cell(n(:,1:3),1);
v = num2cell(v(:,1:3),1);
% n = num2cell( snormalize(rand(3)),1 );
% v = num2cell( snormalize(rand(2,3)),1 );

X = cell(1,3); N = cell(1,3);
for k=1:3
%     X{k} = kron( m(:,k), n(:,k) );
%     N{k} = null(X{k}');
    N{k} = null( kron( v{k}', n{k}' ) );
end

a = CSimLidar;

xy = null( [N{1} N{2}] );
x = snormalize( xy(1:size(N{1},2),:) ) % Normalize x vectors to keep orthonormal nullspaces
N12 = N{1}*x;
xy = null( [N12 N{3}] );
x = snormalize( xy(1:size(N12,2),:) );
N123 = N12*x;
% (*) Do not admit decoupling from N123
% Divide common nullspace
Na = N123(1:3,:);
Nb = N123(4:6,:);
% Normalize nullspace vectors
% Om = N_1' * N_2;
% GT value: lam such that r = N123 * lam
lam = N123' * r;

% Simultaneously orthogonalize N_1 and N_2
Qa = Na'*Na;
Qb = Nb'*Nb;
B  = Na'*Nb;
% Try to simultaneously diagonalize:
[U,~] = eig( Qa );
Da = U' * Qa * U; 
Db = U' * Qb * U;
B = U' * 0.5*(B+B') * U;
% GT value: rho such that lam = U * rho
rho = U' * lam;

tic
[x y z] = sw_Zhou_alt(diag(Da), diag(Db), diag(B), [B(2,3),B(1,3),B(1,2)] );
toc

% Symbols
Srho = sym('r%d',[1 3])';
% sym( 'Srho', 'real' )
Srho' * Da * Srho
Srho' * Db * Srho
Srho' * B * Srho

BB = [ diag(B) ;
       B(1,2) + B(2,1)
       B(2,3) + B(3,2)
       B(3,1) + B(1,3) ];
NBB = null( BB' );

% Impossible to bi-orthonormalize!?
Q = Na'*Na - Nb'*Nb;
% [u,s,v] = svd( Q );
% D = u'*Q*u; % Diagonal matrix (undefined)
[U,D] = eig( Q );
Da = U'*Qa*U;
Db = U'*Qb*U;
% Get a 3-vector base of equation rho'*D*rho = 0 with equal diagonals
Q1 = Na'*Na;
Q2 = Nb'*Nb;
rho_z = sqrt( -D(1,1)*1/D(end) );
rho{1} = [1 0 rho_z]';
delta{1} = U * rho{1};
rho_z = sqrt( -D(2,2)*1/D(1) );
rho{2} = [0 1 rho_z]';
delta{2} = U * rho{2};
% delta = u * rho