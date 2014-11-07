function R = optimizeRotation_NonWeighted( obj )
% Function to optimize with LM the Rotation wrt given correspondences

R = obj.R;
[ R, err, errNorm, W ] = LM_Man_optim(@LM_Fun,R,...
            'space','SO(3)','weighted',false,...
            'debug',obj.debug_level,...
            'maxIters', obj.maxIters,...
            'minChange', obj.minParamChange,...
            'minErrorChange', obj.minErrorChange);
       
    function [residual, jacobian] = LM_Fun( R )
        residual = obj.FErr_Orthogonality( R );
        jacobian = obj.FJac_Orthogonality( R );
    end
        
end

