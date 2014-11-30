function t = optimizeTranslation_2D_NonWeighted( obj, R )
% Function to optimize with LM the Rotation wrt given correspondences

t = obj.t;

L = obj.cam_L;
q = obj.LRF_Q;
K = obj.K;
LM_Fun = @(t)f_Corner_2D( t, R, L, q, K);

[ t, err, errNorm, W ] = LM_Man_optim(LM_Fun,t,...
            'space','Rn','weighted',false,...
            'debug',obj.debug_level,...
            'maxIters', obj.maxIters,...
            'minChange', obj.minParamChange,...
            'minErrorChange', obj.minErrorChange);
         
end