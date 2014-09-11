function t = optimizeTranslation_2D_NonWeighted( obj )
% Function to optimize with LM the Rotation wrt given correspondences

t0 = obj.t0;

R = obj.R;
L = obj.cam_L;
q = obj.LRF_Q;
K = obj.K;
LM_Fun = @(t)f_Corner_2D( t, R, L, q, K);

[ t, err, errNorm, W ] = LM_Man_optim(LM_Fun,t0,...
            'space','Rn','weighted',false,'debug',obj.debug_level, 'maxIters', 200);
        
obj.t = t;     
end