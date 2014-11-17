function [ X, Err ] = Newton_solve( Fun, X0, varargin )
% [ X, Err ] = Newton_solve( Fun, X0, 'option', opt_val, ... )
% Input:
% Fun - function whose outputs are [F, JF]
% X0 - initial values for SO(3), SE(3) or Rn variable
%
% List of options:
% {'debug'; 'space' ; 'maxIters' ; 'minChange' ; 'minErrorChange'};
%
% Space can be (string):
% Manifold - a manifold object
% SO(3) - rotation in 3D
% SE(3) - Euclidean transformation in 3D
% Rn - Euclidean space
%
% Output:
% X - SO(3), SE(3) or R3 transformation
% errNorm - LS error

allfields = {'debug', 1
    'space', 'Rn'
    'maxIters', 15
    'minChange', eps
    'minErrorChange', 1e-14
    };
opt = setVarargin(allfields, varargin{:});
% Algorithm options
maxIters = opt.maxIters; minChange = opt.minChange; minErrorChange = opt.minErrorChange;

iters = 0;
X  = X0;
while iters < maxIters
    iters = iters+1;
    
    [F,JF] = Fun(X);
    step( );
    
    Err = X-X0;
    dif = norm(Err, 'fro');
    if opt.debug
        fprintf('Fro norm = %e\n', dif )
    end
    if dif > minChange
        X0 = X;
    else
        break
    end
end
if opt.debug
    fprintf('Finished in %d iterations with %e final change in Frobenius norm\n', iters,dif)
    disp('End')
end

% Nested function
    function step( )
        % Computes one step in Newton method
        % The update formula dependends on topological space (manifold)
        epsInc = -JF\F;

        switch opt.space
            case 'Manifold'
                % Use general framework for manifolds
                [~,X] = plus( X0, epsInc );
            case 'SE(3)'
                rotInc = rotation_angle_axis( epsInc(1:3) );
                tInc = epsInc(4:6);
%                 rotInc = rotation_angle_axis( epsInc(4:6) );
%                 tInc = epsInc(1:3);
                X(1:3,1:3) = rotInc * X0(1:3,1:3);
                X(1:3,4) = tInc + X0(1:3,4);
            case 'SO(3)'
                rotInc = rotation_angle_axis( epsInc );
                X = rotInc * X0;
            case 'Rn'
                tInc = epsInc;
                X = tInc + X0;
        end
    end
end
