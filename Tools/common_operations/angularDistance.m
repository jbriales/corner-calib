function d = angularDistance( R1, R2 )
% d = angularDistance( R1, R2 )
% d is angular distance given in [deg]
% 
% Angular Distance. Any rotation in SO(3) can be expressed as a rotation through a given angle θ about some axis. The
% angle can always be chosen such that 0 ≤ θ ≤ π, if necessary by reversing the direction of the axis.
% [Rotation Averaging, 4 Distance Measures on SO(3), Richard Hartley]

if nargin == 1
    R2 = eye(3);
end

chordal_distance = norm( R1-R2, 'fro' );
d = 2 * asind( chordal_distance / (2*sqrt(2)) );