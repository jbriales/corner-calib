% Solve the translation problem calculating the uncertainty
%-----------------------------------------------------------
R  = R0; % Initial estimate employed for rotation optimization
%R  = R_c_s_w;
%R  = R_c_s_nw;

L  = [co.L_P2];
q  = cell2mat([co.q]);  %TODO: different sizes of q
K  = img.K;
A_lh = [co.A_lh]; % TODO: is this covariance the good one?
A_q = [co.A_q];

% for i = 600:610
%     L(:,i)    = L(:,i-5) + rand(3,1) * 0.01;
%     q(:,i)    = q(:,i-5) + rand(2,1) * 0.01;
%     A_q(:,i)  = A_q(:,i-5);
%     A_lh(:,i) = A_lh(:,i-5);
% end

% If noise (Blender's dataset)
%w = skewcoords(R);
%w = w + rand([3 1])*0.01;
%R = fast_skewexp(rand(3)*0.01);
% R = R * fast_skewexp(rand(3,1)*0.001);
% t0 = t0.*rand(3,1)*0.001;
%clc
% Linear solution
% for i = 1:N
%     q_norm(:,i) = (q(:,i)/norm(q(:,i)));
%     L_norm(:,i) = (L(:,i)/norm(L(:,i)));
% end
% b = - dot( L, R(:,1:2) * q_norm, 1 )';
% A = L';
% t_lin2 = A \ b;
    
% Linear solution
b = - dot( L, R(:,1:2) * q, 1 )';
A = L';
t_lin = A \ b;
t0 = zeros(3,1);
% LM optimisation without weighting
F_Lev                               = @(t)f_Corner_2D( t, R, L, q, K);
[ t_c_s_nw, err_nw, errNorm, ~ ]    = LM_Man_optim(F_Lev, t0,'space','Rn','debug',2, 'maxIters', 200);
% LM optimisation with weighting
F_Lev_w                             = @(t)f_Corner_2D_W( t, R, L, q, K, A_lh, A_q);
[ t_c_s_w, err_w, errNorm, ~ ]      = LM_Man_optim(F_Lev_w, t0,'space','Rn','debug',2, 'maxIters', 200,'weighted',true);
% LM optimisation without weighting and ERODE
F_Lev_w                             = @(t)f_Corner_2D( t, R, L, q, K);
[ t_c_s_enw, err_enw, errNorm, ~ ]  = LM_Man_optim(F_Lev, t0,'space','Rn','debug',2, 'maxIters', 200,'erode',true);
% LM optimisation with weighting and ERODE
F_Lev_w                             = @(t)f_Corner_2D_W( t, R, L, q, K, A_lh, A_q);
[ t_c_s_ew, err_ew, errNorm, ~ ]    = LM_Man_optim(F_Lev_w, t0,'space','Rn','debug',2, 'maxIters', 200,'erode',true,'weighted',true);

% % Filter the outliers and refine (does not work)
% k = 1;
% clearvars L_e q_e A_lh_e A_q_e;
% for i = 1:size(err_e,1)
%     
%     if abs(err_e(i)) < 10.0*median(err_e)
%         L_e(:,k)    = L(:,i);
%         q_e(:,k)    = q(:,i);
%         A_lh_e(:,k) = A_lh(:,i);
%         A_q_e(:,k)  = A_q(:,i);
%         k = k+1;
%     end    
%     
% end
% F_Lev_w = @(t)f_Corner_2D_W( t, R, L_e, q_e, K, A_lh_e, A_q_e);
% [ t_e1, err_e1, errNorm, ~ ] = LM_Man_optim(F_Lev_w, t0,'space','Rn','debug',2, 'maxIters', 200,'erode',true);

fprintf('Linear (weighted) solution t:\n');
disp(t_lin');
fprintf('Optimised (non-weighted) solution t:\n');
disp(t_c_s_nw');
fprintf('Optimised (weighted) solution t:\n');
disp(t_c_s_w');
fprintf('Robust (non-weighted) solution t:\n');
disp(t_c_s_enw');
fprintf('Robust (weighted) solution t:\n');
disp(t_c_s_ew');


