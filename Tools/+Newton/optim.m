function X = optim( F, JF, X0, space, debug )
% X = optim( F, JF, X0, space, debug )

opt.space = space; % SO(3), SE(3), R3
if ~exist('debug','var')
    debug = 0;
end

if ~exist('space','var')
    switch size(X0,2)
        case 1
            space = 'R3';
        case 3
            space = 'SO(3)';
        case 4
            space = 'SE(3)';
        otherwise
            error('Incorrect dimension for X0')
    end
end

it = 0;
err_tol = eps;
X  = X0;
while it < 50
    it = it+1;
    
    Newton.step
    
    dif = norm(X-X0, 'fro');
    if debug
        fprintf('Fro norm = %e\n', dif )
    end
    if dif > err_tol
        X0 = X;
    else
        break
    end
end
if debug
    fprintf('Finished in %d iterations with %e final change in Frobenius norm\n', it,dif)
    disp('End')
end
