function stereoCompareBumbGT()

S_L = load(fullfile(pwd,'Stereo_Comp','_L.mat'));
S_R = load(fullfile(pwd,'Stereo_Comp','_R.mat'));

% Bumblebee2 GT values
R_L_R_gt = eye(3);
t_L_R_gt = [0.12 0 0]';

% Comparison for complete error function:
R_L_R = S_L.R_ort * S_R.R_ort';
t_L_R = S_L.t_3D - S_R.t_3D;
fprintf('Angular error (deg) for complete: %f\n',angularDistance(R_L_R,R_L_R_gt));
fprintf('Translation error (cm) for complete: %f\n',norm(t_L_R-t_L_R_gt)*100);

% Comparison for 3D global error function:
R_L_R = S_L.R_global_3D * S_R.R_global_3D';
t_L_R = S_L.t_global_3D - S_R.t_global_3D;
fprintf('Angular error (deg) for global: %f\n',angularDistance(R_L_R,R_L_R_gt));
fprintf('Translation error (cm) for global: %f\n',norm(t_L_R-t_L_R_gt)*100);


end