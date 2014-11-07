function obj_LP2 = auxFunction( obj_xi, K )

Ort = [0 -1 ; 1 0];
l = Manifold.P2.empty(0,3);
for k=1:3
    n = Ort * obj_xi.v(k);
    d = - n' * obj_xi.c;
    l(k) = Manifold.P2( K' * [ n ; d ] ); % Camera->Image line conversion is K^(-T), so inverse is applied
end
obj_LP2 = Manifold.Dyn( l(1), l(2), l(3) );
end