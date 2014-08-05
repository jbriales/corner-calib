function [mu, lam, cam_t_w_s, t_c_w, cam_Ms] = solveCamLidarUnknowns( c, q, R_c_w_k, R_c_s, t_c_s )
% [mu, lam] = solveCamLidarUnknowns( c, q, R_c_w_k, R_c_s, t_c_s )
%   c - corner center in calibrated image
%   q - LIDAR point in k-th world axis
%   R_c_w_k - k-th column of R_c_w (rotation matrix world seen from camera)
%   R_c_s - calibration rotation
%   t_c_s - calibration translation
% 
% [mu, lam, cam_t_w_s, t_c_w] = solveCamLidarUnknowns( c, q, R_c_w_k, R_c_s, t_c_s )
% Returns also:
%   - Translation vector from World to Lidar in Camera coordinates
%   - Translation vector from Camera to World

A = [ R_c_w_k -makehomogeneous(c) ];
b = R_c_s(:,1:2) * q + t_c_s;

x = A \ b;

mu  = x(1);
lam = x(2);

if nargout > 2
    cam_t_w_s = mu * R_c_w_k - R_c_s(:,1:2) * q;
    t_c_w = lam * makehomogeneous( c );
    cam_Ms = [ R_c_s(:,1:2) , cam_t_w_s ];
    
    % After this only do:
    %      R_c_w' * cam_Ms
    %      R_c_w' * cam_t_w_s
end

end
