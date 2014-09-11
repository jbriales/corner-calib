function [ X, Err, errNorm, W ] = LM_Man_optim( Fun, X0, varargin )
% [ X, errNorm ] = LM_Man_optim( Fun, X0, 'options' )
% Input:
% Fun - function whose outputs are [F, JF]
% X0 - initial values for SO(3), SE(3) or Rn variable
%
% List of options:
% {'debug'; 'space'; 'f_Err'; 'f_Jac' ; 'maxIters' ; 'minChange' ; 'minErrorChange'};
%
% Space can be:
% SO(3) - rotation in 3D
% SE(3) - Euclidean transformation in 3D
% Rn - Euclidean space
%
% Output:
% X - SO(3), SE(3) or R3 transformation
% errNorm - LS error
%
% The code has been modified from OpenCV implementation of CvHomographyEstimator::refine
% Optional input arguments
allfields = {'debug', 1
    'space', 'Rn'
    'maxIters', 15
    'minChange', eps
    'minErrorChange', 1e-14
    'weighted', false
    'erode', false
    };
opt = setVarargin(allfields, varargin{:});
% Algorithm options
maxIters = opt.maxIters; minChange = opt.minChange; minErrorChange = opt.minErrorChange;
DONE = 0; STARTED = 1; CALC_J = 2; CHECK_ERR = 3;
Cstate = {'DONE','STARTED','CALC_J','CHECK_ERR'};
lambda = 1e-3; iters = 0; state = STARTED;
% Npoints = length(X);
% -------------------------------------------------------------------------
% Algorithm
% -------------------------------------------------------------------------
prevParam = X0; param = X0;
prevErrNorm = realmax('double');
while(1) % The loop is stopped from inside with internal conditions for convergence
    exit_optim = updateAlt( );
    if ~exit_optim
        X = param;
        if ~exist('W','var')
            W = eye(length(Err));
        end
        % if opt.debug >= 2
        % % Temporal debug
        % figure
        % subplot(211)
        % bar(Err)
        % subplot(212)
        % bar(diag(W))
        % fprintf('Final error: %.3f\n',norm(Err))
        % fprintf('Final weighted error: %.3f\n',errNorm)
        % end
        break
    end
    % WARNING: COULD NEED NEGATIVE SIGN (DEPENDING ON ERROR DEFINITION)
    % See Lev-Marq in Wikipedia
    % step update is positive if err = y_i - Fun(x_i,param)
    % Here Fun(x_i,param) is equivalent to f_Err so Err = 0 - f_Err (values
    % y_i are zero)
    if ~opt.weighted && ~opt.erode
        [Err, J] = Fun( param );
        Err = -Err; % Correct sign (see above)
        M = 1;
        JtJ = J'*J;
        JtErr = J'*Err;
        errNorm = Err'*Err/M;
    else if opt.erode && ~opt.weighted
            [Err, J] = Fun( param );
            Err = -Err; % Correct sign (see above)
            M = length(Err);
            b = 3; b2 = b^2; %Huber function parameter
            JtJ = J'*J;
            F_r = Err.^2;
            sqr = sqrt(1+(F_r/b2));
            errNorm = sum( b2*(sqr+1) )/2;
            d_rho = sum(1./sqr);
            JtErr = d_rho*(J'*Err)/M;
        else if opt.erode && opt.weighted
                [Err, J, W] = Fun( param );
                Err = -Err; % Correct sign (see above)
                M = length(Err);
                b = 3; b2 = b^2; %Huber function parameter
                JtJ = J'*W*J;
                F_r = W*Err.^2;
                sqr = sqrt(1+(F_r/b2));
                errNorm = sum( b2*(sqr+1) )/2;
                d_rho = sum(1./sqr);
                JtErr = d_rho*(J'*W*Err)/M;
            else
                [Err, J, W] = Fun( param );
                Err = -Err; % Correct sign (see above)
                M = 1;
                JtJ = J'*W*J;
                JtErr = J'*W*Err;
                errNorm = Err'*W*Err/M;
            end
        end
    end
    if opt.debug>=1
        fprintf('It = %d\t|Err|=%e\tState is %s\tlog10(lambda)=%d\n',iters,errNorm,Cstate{state+1},log10(lambda))
    end
end
if opt.debug>=1
    if exit_cond(1)
        fprintf('maxIters reached: iters = %d\tmaxIters = %d\n',iters,maxIters)
    end
    if exit_cond(2)
        fprintf('paramChange below minChange: paramChange = %e\tminChange = %e\n',paramChange,minChange)
    end
    if exit_cond(3)
        fprintf('errorChange below minErrorChange: errorChange = %e\tminErrorChange = %e\n',errorChange,minErrorChange)
    end
    fprintf('\n')
end
%% Nested functions
    function exit_optim = updateAlt( )
        % Update the state
        if( state == DONE )
            exit_optim = false;
            return
        end
        if( state == STARTED )
            state = CALC_J;
            exit_optim = true;
            return
        end
        if( state == CALC_J )
            prevParam = param;
            step()
            prevErrNorm = errNorm;
            state = CHECK_ERR;
            exit_optim = true;
            return
        end
        if ( state ~= CHECK_ERR )
            warning('Bug in state values');
        end
        if( errNorm > prevErrNorm )
            lambda = lambda * 10;
            if( log10(lambda) <= 16 )
                step()
                state = CHECK_ERR;
                exit_optim = true;
                return
            end
        end
        if log10(lambda) >= -16
            lambda = lambda / 10;
        end
        iters = iters + 1;
        % TODO: Measure paramChange in manifold
        paramChange = norm( param - prevParam, 2 );
        errorChange = prevErrNorm - errNorm;
        if errorChange < 0
            warning('Error change in last LM iteration is negative (increased error)');
            warning('Error function could not be decreased: Check error function is correct');
            exit_optim = false;
            exit_cond = false(1,3);
            return
        end
        exit_cond = [iters > maxIters , paramChange < minChange , errorChange < minErrorChange];
        if( any( exit_cond ) )
            state = DONE;
            exit_optim = false;
            return
        end
        prevErrNorm = errNorm;
        state = CALC_J;
        exit_optim = true;
    end
    function step( )
        % Computes one step in Levenberg-Marquardt
        % The update formula dependends on topological space (manifold)
        JtJN = JtJ + lambda * diag(diag(JtJ));
        epsInc = JtJN \ JtErr;
        % WARNING: COULD NEED NEGATIVE SIGN (DEPENDING ON ERROR DEFINITION)
        switch opt.space
            case 'SE(3)'
                % rotInc = rotation_angle_axis( epsInc(1:3) );
                % tInc = epsInc(4:6);
                rotInc = rotation_angle_axis( epsInc(4:6) );
                tInc = epsInc(1:3);
                param(1:3,1:3) = rotInc * prevParam(1:3,1:3);
                param(1:3,4) = tInc + prevParam(1:3,4);
            case 'SO(3)'
                rotInc = rotation_angle_axis( epsInc );
                param = rotInc * prevParam;
            case 'Rn'
                tInc = epsInc;
                param = tInc + prevParam;
        end
    end
end
function opt = setVarargin( allfields, varargin )
% opt = setVarargin( allfields, varargin )
% Input:
% - Nx2 cell of optional parameters
% - varargin cell with pairs of option-value parameters
opt = struct();
nVarargs = length(varargin);
if mod(nVarargs,2)~=0
    error('Optional arguments should be pairs: "Property name", Value')
end
for k = 1:2:nVarargs
    opt.(varargin{k}) = varargin{k+1};
end
% Default values in options
for k = 1:size(allfields,1)
    if ~isfield(opt,allfields{k,1})
        opt.(allfields{k,1}) = allfields{k,2};
    end
end
end