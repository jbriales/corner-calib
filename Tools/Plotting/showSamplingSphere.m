function showSamplingSphere(N,L)

vecL = zeros(3, size(N,2));
for i=1:size(N,2) % For representing measures in 3D space
    n = N(:,i);
%     J_n_y = Householder(n) * [eye(2) , [0 0]']';
    J_n_y = Householder(n);
    J_n_y = J_n_y(:,1:2);
    vecL(:,i) = J_n_y * L(:,i);
end
figure('Name','Sampling map'), hold on
scatter3( N(1,:), N(2,:), N(3,:), 'r' )
quiver3( N(1,:), N(2,:), N(3,:), vecL(1,:), vecL(2,:), vecL(3,:), 0.1, 'b' )
axis equal, rotate3d on

end

function H = Householder( x, ei )
% H  = Householder( x, ei )
%   H is the Householder matrix associated to x

if ~exist('ei','var')
    ei = 3;
end
n = length( x );
I = eye( n );

v = x - sign( x(1) ) * norm(x) * I(:,ei);

H = I - 2 * v * v' / (v'*v);

end