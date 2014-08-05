% Basis for simulation main

% Generate random samples
% Should give:
% For image - img_params (p,ang1,ang2,ang3) in pixel values and its
% covariance
% For LIDAR - scan pts and cell array with 3 groups of inliers
generate_random_poses
    
% Functions to apply to images
[NL, c, ~, L_P2] = imgParams2calibratedData( img_params, img.K );
[R_c_w, A_N, A_eps] = getWorldNormals( R0, NL, c, A_co );

% Functions to apply to LIDAR
[l,A_l,A_lh,p,A_p,q,A_q, lin,seg] = computeScanCorner( scan, {scantrack.all_inliers}, debug );