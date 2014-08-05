function gt = loadGT( gt, path, im_frame )

[~,v] = collectPoses_Blender(path);
% Copied from Blender
gt.Rx = v(5,im_frame);
gt.Ry = v(6,im_frame);
gt.Rz = v(7,im_frame);

R_w_b = RotationZ( deg2rad( gt.Rz ) ) * RotationY( deg2rad( gt.Ry ) ) * RotationX( deg2rad( gt.Rx ) );
R_b_c = diag([1 -1 -1]);
R_w_c = R_w_b * R_b_c;
gt.R  = R_w_c'; % gt.R is R_c_w: World pose seen from Camera

gt.R_w_c = R_w_c;
gt.R_c_w = R_w_c';
gt.R_w_s = gt.R_w_c * gt.R_c_s;

gt.tx = v(2,im_frame);
gt.ty = v(3,im_frame);
gt.tz = v(4,im_frame);

gt.t_w_c = [gt.tx gt.ty gt.tz]';
gt.t_c_w = -gt.R_w_c' * gt.t_w_c; % gt.t_c_w is World pose seen from Camera
gt.t_w_s = gt.t_w_c + gt.R_w_c * gt.t_c_s; % negative ? No, positive

end