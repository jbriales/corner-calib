function [obj_Nbp, obj_LP2] = computeBackprojectedNormals( obj_xi, K )
% [obj_Nbp, obj_LP2] = getBackprojectedNormals( obj_xi, K )
% Extract useful information from image xi parameter
% Input:
%   - xi object (TO image parameters)
%   - K camera intrinsic calibration
% Output:
%   - obj_Nbp object with back-projected planes normals
%   - obj_LP2 object with homogeneous lines in image

% Compute homogeneous lines (in camera projective space)
Ort = [0 -1 ; 1 0];
l = Manifold.P2.empty(0,3);
for k=1:3
    n = Ort * obj_xi.v(k);
    d = - n' * obj_xi.c;
%     l{k} = K' * [ n ; d ];
    l(k) = Manifold.P2( K' * [ n ; d ] ); % Camera->Image line conversion is K^(-T), so inverse is applied
end
J_project_P2 = blkdiag( l(1).DLie, l(2).DLie, l(3).DLie );
J_L_xi = J_project_P2 * Jacobian_L_xi ( obj_xi, K );
A_L = J_L_xi * obj_xi.A_X * J_L_xi'; % Covariance of 3 P2 lines in camera frame
obj_LP2 = Manifold.Dyn( l(1), l(2), l(3) );
obj_LP2.setRepresentationCov( A_L );

% Compute back-projected planes normals (final feature)
J_normalize = cell(1,3);
N = Manifold.S2.empty(0,3);
for k=1:3
    N(k) = Manifold.S2( l(k).X / norm(l(k).X) );
    J_normalize{k} = N(k).DLie * 1/norm(l(k).X)^3 * (eye(3) - l(k).X*l(k).X');
end
J_normalize = blkdiag( J_normalize{:} );
A_N = J_normalize * A_L * J_normalize';

obj_Nbp = Manifold.Dyn( N(1), N(2), N(3) );
obj_Nbp.setRepresentationCov( A_N );

end

function [ J_L_xi ] = Jacobian_L_xi ( obj_xi, K )
% Jacobian of 3 camera P2 lines wrt xi object (c,v1,v2,v3)
Ort = [0 -1 ; 1 0];
c  = obj_xi.c;
v1 = Ort * obj_xi.v(1).X;
v2 = Ort * obj_xi.v(2).X;
v3 = Ort * obj_xi.v(3).X;

Jc = @(v) K' * [ zeros(2,2) ; -(Ort*v)' ];
Jv = cell(1,3);
for i=1:3
    Jv{i} = K' * [ eye(2) ; -c' ] * Ort * obj_xi.v(1).DLie; % Projected to S1
end
J_L_xi = [ [Jc(v1); Jc(v2); Jc(v3)] , blkdiag(Jv{:}) ];

end