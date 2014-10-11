% All scripts

% Global variables
global WITH_MONTECARLO %#ok<NUSED>

% TODO: Create two auxiliary functions for tracking and storing
% metadata in images and scans before optimization

% Previous assertions:
clear
% close all
dbstop if error % Set debug mode on if error occurs
kbhit('stop')
kbhit('init') % Track keyboard hits

% To position figures centered and visible by default
set(0,'defaultfigureposition', [642   503   560   420])

main_type = readConfigFile( fullfile(pwd,'main.ini'),'','main_type' );
userOpts = readConfigFile( fullfile(pwd,'main.ini'), main_type );
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

            % Temporal
            WITHGT = false;
            
%% FOR LOOP BEGIN
for nobs=1:length(scans)
    fprintf('nframe = %d\n',nobs)
    fprintf('-------------\n')
    figure( hWin )
    
    if upper(kbhit) == 'D' % To stop after pressing 'L' key
        disp('D was hit: Debug activation')
        keyboard
    end

    img = imgs(nobs);
    
    if WITHGT
        gt = loadGT( gt, path, nobs );
        co(nobs).gt = gt;
    end
    
    %% Trihedron data from image
    % TO image is named as xi, composed of:
    %   c - vertex image and 3 lines, R2
    %   v(i) - 3 directions in image, S1
    % Load image: Remove if not necessary to improve speed
    Cam.setFrame( frames(nobs) );
    obj_xi = Cam.computeXi;
    if strcmp( win_visibility, 'on' )
        Cam.visualize;
    end
    
    [obj_Nbp, obj_LP2] = computeBackprojectedNormals( obj_xi, Cam.K );
    
    % Compute trihedron planes normals from xi parameter
    tic
    obj_Rtri = computeTrihedronNormals( obj_xi, Cam.K, obj_Nbp );
    fprintf('R_c_w optimization TIME: %f\n',toc)
    if WITHGT
        disp('Compare computed R_c_w to GT')
        disp(obj_Rtri.X), disp(gt.R_c_w)
        fprintf('Angular distance with GT rotation: %f\n', angularDistance(obj_Rtri.X, gt.R_c_w))
        if angularDistance(obj_Rtri.X, gt.R_c_w) > 3
            CHECK_IMAGE = true;
        end
    end
    
    % Check displacement from previous frame
    if exist('R_c_w0','var')
        fprintf('Angular distance from image %d to %d: %e\n', nobs-1, nobs, angularDistance(R_c_tri0,obj_Rtri.X))
    end
    R_c_tri0 = obj_Rtri.X; % Store previous for comparison
    
    debug_im_reprojection = false;
    if debug_im_reprojection
        Cam.visualizeTrihedronReprojection;
    end
    
    %% Go from LIDAR information to Corner Observation
    % Set algorithm input
    LRF.setFrame( scans(nobs) );
    [v,A_v,q,A_q] = LRF.compute;
    if strcmp( win_visibility, 'on' )
        LRF.visualize;
    end
    
    % User reallocation of LIDAR segments
    if upper(kbhit) == 'L' % To stop after pressing 'L' key
        disp('L was hit: Manual reallocation of segments')
        keyboard
        
        if ~exist('thereis_line','var') || all(thereis_line)
            selection_mask = input('Input segment selection mask [x y z]: ');
        else
            selection_mask = ~thereis_line;
        end
        
        subplot( hLidar ), set( gcf, 'Visible', win_visibility )
        scantrack_aux = manualSetScanlines( scan.xy, selection_mask );
        for k=1:3
            if selection_mask(k) % Substitutes only user selected segments
                scantrack(k) = scantrack_aux(k);
            end
        end
    end
    
%     debug = 0;
%     [v,A_v,~,~,q,A_q, lin,seg] = ...
%         computeScanTO( LRF.xy, LRF.sd,...
%         {LRF.meta.scantrack.all_inliers}, debug );
%     [v,A_v,~,~,q,A_q, lin,seg] = computeScanTO( scan.xy, scan_sigma, {scantrack.all_inliers}, debug );
%     [l,A_l,A_lh,p,A_p,q,A_q, lin,seg] = computeScanCorner( scan, scan_sigma, {scantrack.all_inliers}, debug );
%     [scantrack.lin] = deal( lin{:} );
%     [scantrack.seg] = deal( seg{:} );    
        
    %% Store data for final optimization
    tic
    co = CTrihedronObservation( obj_Rtri, obj_LP2, obj_Nbp,...
                v, A_v, [], [], q, A_q );
	triOptim.stackObservation( co );
    fprintf('STORAGE TIME: %f\n',toc)
    
end
save( fullfile( path, 'cache',...
      strcat('preoptimization',stereoLabel,hokuyoLabel) ),...
      'triOptim' )
keyboard
figure('Name','Sampling map'), hold on
triOptim.showSamplingSphere;

%% Final optimization with triOptim object
% WITHRANSAC = false;
% WITHRANSAC = true;

if 0 % Code to filter obs for R and t from different datasets
    mask_data1 = kron([frames.ts]<1412150000, ones(1,3));
    triOptim.mask_RANSAC_t_outliers = mask_data1;
    mask_data2 = kron([frames.ts]>1412150000, ones(1,3));
    triOptim.mask_RANSAC_R_outliers = mask_data2;
end

if WITHRANSAC
    triOptim.filterRotationRANSAC;
    triOptim.showSamplingSphere('r');
else % Give initial estimate
%     triOptim.setInitialRotation( gt.R_c_s );
    triOptim.setInitialRotation( [ 0 -1 0;
                                   0 0 -1;
                                   +1 0 0] );
end
triOptim.disp_N_R_inliers;
R_ort = triOptim.optimizeRotation_Weighted;

counter = 0; dist = +inf;
while dist > 1e-2 && counter < 10
    triOptim.resetRotationRANSAC;
    triOptim.filterRotation_R0( R_ort );
    triOptim.disp_N_R_inliers;
    disp(R_ort)
    R_ort_0 = R_ort;
    R_ort = triOptim.optimizeRotation_Weighted;
    dist = angularDistance( R_ort, R_ort_0 );
    
    counter = counter + 1;
end
if counter==10
    warning('Counter %d reached in rotation',counter);
end

if WITHRANSAC
    triOptim.filterTranslationRANSAC( R_ort );
else
%     triOptim.setInitialTranslation( gt.t_c_s );
    triOptim.setInitialTranslation( [ 0.4 0 0.15 ]' );
end
triOptim.disp_N_t_inliers;
% triOptim.setInitialTranslation( Rig.t_c_s + 0.05*randn(3,1) );
% t_3D_nw = triOptim.optimizeTranslation_3D_NonWeighted( R_c_s_w );
t_3D = triOptim.optimizeTranslation_3D_Weighted( R_ort );
% t_2D_nw = triOptim.optimizeTranslation_2D_NonWeighted( R_c_s_w );
% t_2D_w = triOptim.optimizeTranslation_2D_Weighted( R_c_s_w );

counter = 0; dist = +inf;
while dist > 1e-3 && counter < 10
    triOptim.resetTranslationRANSAC;
    triOptim.filterTranslation_t0( R_ort, t_3D );
    triOptim.disp_N_t_inliers;
    disp(t_3D)
    t_3D_0 = t_3D;
    t_3D = triOptim.optimizeTranslation_3D_Weighted( R_ort );
    dist = norm( t_3D-t_3D_0 );
    
    counter = counter + 1;
end
if counter==10
    warning('Counter %d reached in translation',counter);
end

[R_global, t_global] = triOptim.optimizeGlobal_Ort_3D( R_ort, t_3D );
[R_global_3D, t_global_3D] = triOptim.optimizeGlobal_3D( R_ort, t_3D );

% Code to reproject LRF to img
T = [ 0 -1 0 ; 0 0 -1 ; 1 0 0 ];
R_m = RotationY(deg2rad(-20)) * T;
% t_m = -R_m * [0 0.54 0]';
t_m = [ 0.5 -0.05 0.20 ]';
Rt_m = [ R_m t_m ];

disp( [R_ort t_3D] )
% if exist('gt','var')
%    disp( [gt.R_c_s gt.t_c_s] )
% end
disp( Rt_m )
disp( [R_global_3D, t_global_3D] )

if stereoStore % Save stereo results for later comparison and study
    filename = fullfile(pwd,'Stereo_Comp',stereoLabel);
    save(filename, 'R_ort', 't_3D', 'R_global_3D', 't_global_3D');
end

if 0
    keyboard
    im_gt = imread(fullfile(pwd,'GT',strcat('CAM',stereoLabel,'_2.jpg')));
%     xy_gt = load(fullfile(pwd,'GT','XY1.mat'),'XY');
%     xy_gt = xy_gt.XY;
    meanLRF % Extracts GT data
    K = [	248.321289 0 319.416809 ;
	0 248.321289 249.839676 ;
	0 0 1 ];
    reprojectScan2Img( im_gt, xy_gt, K, [R_m t_m] );
    reprojectScan2Img( im_gt, xy_gt, K, [R_ort t_3D] );
    reprojectScan2Img( im_gt, xy_gt, K, [R_global_3D t_global_3D] );
%     [gt.R_c_s gt.t_c_s]
end
if 0 % Check GT images
    im_gt = imread('/media/storage/Research/Corner_Blender/Corner_Comercio_GT/Im1.png');
    xy_gt = loadBlenderPCD('/media/storage/Research/Corner_Blender/Corner_Comercio_GT/Scan1.pcd');
    reprojectScan2Img( im_gt, xy_gt, Cam.K, [R_m t_m] );
    reprojectScan2Img( im_gt, xy_gt, Cam.K, [R_ort t_3D], false );
    reprojectScan2Img( im_gt, xy_gt, Cam.K, [R_global_3D t_global_3D], false );
end
return