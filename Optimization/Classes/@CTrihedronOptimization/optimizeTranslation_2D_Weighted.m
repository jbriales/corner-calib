function t = optimizeTranslation_2D_Weighted( obj )
% Function to optimize with LM the Rotation wrt given correspondences

t0 = obj.t0;

R = obj.R;
L = obj.cam_L;
q = obj.LRF_Q;
K = obj.K;
% TODO: Read covariances filtering with RANSAC
LM_Fun = @(t)f_Corner_2D_W( t, R, L, q, K, A_lh, A_q);

[ t, err, errNorm, W ] = LM_Man_optim(LM_Fun,t0,...
            'space','Rn','weighted',false,'debug',obj.debug_level, 'maxIters', 200);
        
obj.t = t;     
end