function gt = loadCalibrationGT()
% Set blender groundtruth

% Get real transformation (C to S in standard SR's) (Grountruth)
gt.R_cb_c = diag([1 -1 -1]);
gt.R_sb_s = [ 0 -1 0 ;
           0  0 1 ;
          -1  0 0 ];

R_cb_sb = RotationX(deg2rad(-7)); % Rotation in Blender angles
t_cb_sb = [0 0.5 0]'; % Translation in Blender frame

R_c_s = gt.R_cb_c' * R_cb_sb * gt.R_sb_s; % Rotation matrix in standard SR
t_c_s = gt.R_cb_c' * t_cb_sb;
Rt_c_s = [R_c_s t_c_s];

gt.R_c_s = R_c_s;
gt.t_c_s = t_c_s;

end