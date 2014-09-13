function t = optimizeTranslation_3D_Weighted( obj, R )
% Function to optimize with LM the translation wrt given correspondences
% and Rotation previous estimate

t0 = obj.t0;
[ t, err, errNorm, W ] = LM_Man_optim(@LM_Fun,t0,...
            'space','Rn','weighted',true,...
            'debug',obj.debug_level,...
            'maxIters', obj.maxIters,...
            'minChange', obj.minParamChange,...
            'minErrorChange', obj.minErrorChange);
       
    function [residual, jacobian, weights] = LM_Fun( t )
        residual = obj.FErr_3D_PlaneDistance( R, t );
        jacobian = obj.FJac_3D_PlaneDistance( R, t );
        weights  = obj.FWeights_3D_PlaneDistance( R, t );
    end
        
end