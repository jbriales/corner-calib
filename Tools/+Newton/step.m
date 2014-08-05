% function step( )
% Computes one step in Newton method
% global param prevParam JtJ JtErr lambda

epsInc = -JF(X)\F(X);
% WARNING: COULD NEED NEGATIVE SIGN (DEPENDING ON ERROR DEFINITION)
% See notes in LevMarq.optim
switch opt.space
    case 'SE(3)'
        rotInc = rotation_angle_axis( epsInc(1:3) );
        tInc   = epsInc(4:6);
        
        X(1:3,1:3) = rotInc * X0(1:3,1:3);
        X(1:3,4)   = tInc   + X0(1:3,4);
    case 'SO(3)'
        rotInc = rotation_angle_axis( epsInc );
        X = rotInc * X0;
    case 'R3'
        tInc = epsInc;
        X = tInc + X0;
end