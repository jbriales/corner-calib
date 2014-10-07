% All scripts

% Global variables
global WITH_MONTECARLO %#ok<NUSED>

% TODO: Create two auxiliary functions for tracking and storing
% metadata in images and scans before optimization

% Previous assertions:
% close all

% To position figures centered and visible by default
set(0,'defaultfigureposition', [642   503   560   420])

userOpts = readConfigFile( fullfile(pwd,'main.ini'), '[GT_IR]' );
extractStructFields( userOpts );

all_dataset_data = loadDataset( userOpts );
extractStructFields( all_dataset_data );

% debug = 0;
% [img_params, A_img_params, imgtrack, checkImage] = ...
%     corner_calib(imgtrack, rgb2gray(imgs(1).I), debug);

CompleteCO = true;
repr_planes = [];
scan_points = [];
CHECK_IMAGE = false;

% Create objects for storage and optimization
optim_config_file = fullfile( pwd, 'optim_config.ini' );
optimOpts = readConfigFile( optim_config_file );
extractStructFields( optimOpts );
clear optimOpts
triOptim = CTrihedronOptimization( Cam.K,...
    RANSAC_Rotation_threshold,...
    RANSAC_Translation_threshold,...
    debug_level, maxIters,...
    minParamChange, minErrorChange);


% For Left
if strcmp(stereoLabel,'_L')
    R_ort_0 = RotationY(deg2rad(-1))* [ -0.3434   -0.9381    0.0443
              -0.0635   -0.0238   -0.9977
               0.9370   -0.3454   -0.0514 ];
    t_3D_0 = [     0.5128    -0.0356     0.1715]';
% For Right
elseif strcmp(stereoLabel,'_R')
    R_ort_0 = [ -0.3481   -0.9365    0.0428
        -0.0589   -0.0237   -0.9980
        0.9356   -0.3499   -0.0469 ];
    t_3D_0 = [    0.3956    -0.0388     0.1739]';
end
% if strcmp(stereoLabel,'_L')
%     R_ort = [ -0.3434   -0.9381    0.0443
%               -0.0635   -0.0238   -0.9977
%                0.9370   -0.3454   -0.0514 ];
%     t_3D = [     0.5128    -0.0356     0.1715]';
% % For Right
% elseif strcmp(stereoLabel,'_R')
%     R_ort = [ -0.3481   -0.9365    0.0428
%         -0.0589   -0.0237   -0.9980
%         0.9356   -0.3499   -0.0469 ];
%     t_3D = [    0.3956    -0.0388     0.1739]';
% end

% Happy idea
R_ort = [  -0.2825   -0.9588    0.0307
           -0.0559   -0.0155   -0.9983
            0.9576   -0.2837   -0.0492 ];
t_3D = [ 0.3356 ; -0.0588 ; 0.0739 ];

%% FOR LOOP BEGIN
% close all
for nobs=1:length(scans)
    fprintf('nframe = %d\n',nobs)
    fprintf('-------------\n')
    figure( hWin )
    
    if upper(kbhit) == 'D' % To stop after pressing 'L' key
        disp('D was hit: Debug activation')
        keyboard
    end
    
    Cam.hFig = figure('Name','Reprojection to camera');
    Cam.frame = frames(nobs);
    Cam.loadImage;
    Cam.showImage;
%     if strcmp( win_visibility, 'on' )
%         Cam.visualize;
%     end
    
    LRF.setFrame( scans(nobs) );
%     LRF.xy = scans(nobs).xy;
    
%     reprojectScan2Img( im_gt, xy_gt, K, [R_m t_m] );
    hold on
%     reprojectScan2Img( Cam.img, LRF.xy, Cam.K, [R_ort t_3D], false );
   reprojectScan2ImgCol( Cam.img, LRF.xy(:,LRF.meta.scantrack(1).inliers{1}),...
                       Cam.K, [R_ort t_3D], 'w' );
   reprojectScan2ImgCol( Cam.img, LRF.xy(:,LRF.meta.scantrack(1).inliers{2}),...
       Cam.K, [R_ort t_3D], 'w' );
   reprojectScan2ImgCol( Cam.img, LRF.xy(:,LRF.meta.scantrack(1).inliers{3}),...
       Cam.K, [R_ort t_3D], 'w' );
   
   % Test
   if 0
       pts = [ LRF.xy(:,LRF.meta.scantrack(1).inliers{1}) ,...
           LRF.xy(:,LRF.meta.scantrack(1).inliers{2}) ,...
           LRF.xy(:,LRF.meta.scantrack(1).inliers{3}) ];
       Npts = size(pts,2);
%        delete(hPts);
       t_3D = t_3D_0 + [ -0.06 ; -0.02 ; -0.10 ];
       R_ort = RotationY(deg2rad(4)) * RotationZ(deg2rad(-0.5)) * R_ort_0;
       pts2D = makeinhomogeneous( Cam.K*( R_ort(:,1:2)*pts + repmat(t_3D,1,Npts) ) );
       hPts = scatter( pts2D(1,:), pts2D(2,:), [], 'w' );
   end

%     reprojectScan2Img( im_gt, xy_gt, K, [R_global_3D t_global_3D] );

%     pause( 0.1 )
end