for k=1:3
    % Create trihedron lines directions
    v{k} = snormalize(rand(3,1));
end
v_ = [v{:}];
   
for k=1:3
    % Compute normal directions
    ij = setdiff(1:3,k); i = ij(1); j = ij(2);
    n{k} = snormalize( cross(v{i},v{j}) );
end
n_ = [n{:}];

figure, hold on
Os = zeros(1,3);
hv = quiver3( Os, Os, Os, v_(1,:), v_(2,:), v_(3,:), 'Color','r' );
hn = quiver3( Os, Os, Os, n_(1,:), n_(2,:), n_(3,:), 'Color','b' );
axis equal, rotate3d on

% Compute scalar products (interior angle)
for k=1:3
    ij = setdiff(1:3,k); i = ij(1); j = ij(2);
    disp( v{i}'*v{j} )
    disp( n{i}'*n{j} )
end

% Compute symbolic expression for angles between trihedron lines wrt 3
% parameters
syms a b c
canvec(1)' * RotationZ( a ) * canvec(2)
canvec(2)' * (RotationZ(c-a) * RotationY(b)) * canvec(3)
canvec(3)' * (RotationY(b).' * RotationZ(c).') * canvec(1)
