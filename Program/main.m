% All scripts

% Global variables
global WITH_MONTECARLO %#ok<NUSED>

% TODO: Create two auxiliary functions for tracking and storing
% metadata in images and scans before optimization

% Previous assertions:
clear
close all
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
triOptim = CTrihedronOptimization( imgs(1).K,...
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
    Cam.visualize;
    
    [obj_Nbp, obj_LP2] = computeBackprojectedNormals( obj_xi, Cam.K );
    
    % Compute trihedron planes normals from xi parameter
    tic
    obj_Rtri = computeTrihedronNormals( obj_xi, Cam.K, obj_Nbp );
%     [R_c_w, A_R_c_w, A_eps_c_w] = getWorldNormals( R0, NL, c, A_co );
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
    scan = scans(nobs);
    if isempty(scan.xy)
        switch typeOfSource
            case 'Blender'
                scan.xy = loadBlenderPCD( scan.path );
            case 'Rawlog'
                % TODO: Not necessary?
        end
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
    
    debug = 0;
    metafile = scan.metafile;
    if exist(metafile,'file')
        % TODO: What to do with lost var?
        load( metafile, '-mat', 'scantrack', 'inPts', 'lost' );
    else
        tic
        [scantrack, inPts, lost] = updateSegments( scan, scantrack, debug );
        save( metafile, '-mat', 'scantrack', 'inPts', 'lost' );
    end
    clear metafile
    [v,A_v,~,~,q,A_q, lin,seg] = computeScanTO( scan.xy, scan_sigma, {scantrack.all_inliers}, debug );
%     [l,A_l,A_lh,p,A_p,q,A_q, lin,seg] = computeScanCorner( scan, scan_sigma, {scantrack.all_inliers}, debug );
    [scantrack.lin] = deal( lin{:} );
    [scantrack.seg] = deal( seg{:} );
    
    fprintf('SCAN2CO TIME: %f\n',toc)
    thereis_line   = cellfun(@(x)~isempty(x), v);
    thereis_corner = cellfun(@(x)~isempty(x), q);
    if all(thereis_corner) && all(thereis_line)
        CompleteCO = true;
    else
        CompleteCO = false;
        if 0 % Automatic reallocation
            % Estimate wrt to thereis_corner data
            [~, ~, ~, ~, cam_Ms] = solveCamLidarUnknowns( c, q{k}, R_c_w(:,k), R_c_s, t_c_s );
            Ms = R_c_w' * cam_Ms;
        end
    end
    
    if any( lost ) % Stop in iteration with loss of complete CO
        lab = 'XYZ';
        warning('LIDAR tracking has lost segment %c', lab(lost))
        keyboard
    end
    
    %% Update LIDAR visualization
    % TODO: Put as option in a new class CLRF
    tic
    subplot(hLidar), set( gcf, 'Visible', win_visibility )
    ax = axis;
    cla, hold on
    hold on, title('Current scan segmentation')
    plot( scan.xy(1,:), scan.xy(2,:), '.k' ), axis equal
    axis( ax );
    plotLIDARframe( );
    col = 'rgb';
    for k=1:3
        if thereis_line(k)
            plot( inPts{k}(1,:), inPts{k}(2,:), [col(k),'.'] );
            plotHomLineWin( scantrack(k).lin, col(k) );
        end
    end
    fprintf('LIDAR PLOT TIME: %f\n',toc)
        
    %% Store data for final optimization
    tic
    co = CTrihedronObservation( obj_Rtri, obj_LP2, obj_Nbp,...
                v, A_v, [], [], q, A_q );
	triOptim.stackObservation( co );
    if 0 % Deprecated, use new classes
        lab = [1 2 3];
        co(nobs).complete = CompleteCO;
        co(nobs).thereis_line = thereis_line; % Line - rotation correspondence
        co(nobs).lab_line = lab(thereis_line);
        co(nobs).thereis_corner = thereis_corner; % Corner - translation correspondence
        co(nobs).lab_corner = lab(thereis_corner);
        
        % For rotation optimization
        co(nobs).R_c_w   = R_c_w; % Normal vectors to world planes (in Camera SR)
        co(nobs).A_R_c_w = A_R_c_w;
        co(nobs).l       = l; % Direction vector of LIDAR segments (in LIDAR SR)
        co(nobs).A_l     = A_l;
        
        % For translation optimization
        co(nobs).N_repr  = N_repr; % Normal vectors to reprojection planes (in Camera SR)
        co(nobs).L_P2    = L_P2;
        co(nobs).A_lh    = A_lh;
        co(nobs).q       = q; % 2D (XY) Corner points (in LIDAR SR)
        co(nobs).A_q     = A_q;
        
        % For final checking and debug
        check(nobs).cam_frame  = nobs;
        check(nobs).scan_frame = nobs;
        check(nobs).scan_inPts = inPts;
        check(nobs).scan_track = scantrack;
        check(nobs).img_track  = imgtrack;
    end
    fprintf('STORAGE TIME: %f\n',toc)
    
end
save( fullfile( path, 'cache',...
      strcat('preoptimization',stereoLabel,hokuyoLabel) ),...
      'triOptim' )
keyboard

%% Final optimization with triOptim object
WITHRANSAC = false;
if WITHRANSAC
    triOptim.filterRotationRANSAC;
else % Give initial estimate
    triOptim.setInitialRotation( gt.R_c_s );
end
triOptim.disp_N_R_inliers;
R_c_s_w  = triOptim.optimizeRotation_Weighted;

if WITHRANSAC
    triOptim.filterTranslationRANSAC( R_c_s_w );
else
    triOptim.setInitialTranslation( gt.t_c_s );
end
triOptim.disp_N_t_inliers;
% triOptim.setInitialTranslation( Rig.t_c_s + 0.05*randn(3,1) );
% t_3D_nw = triOptim.optimizeTranslation_3D_NonWeighted( R_c_s_w );
t_3D_w  = triOptim.optimizeTranslation_3D_Weighted( R_c_s_w );
% t_2D_nw = triOptim.optimizeTranslation_2D_NonWeighted( R_c_s_w );
% t_2D_w = triOptim.optimizeTranslation_2D_Weighted( R_c_s_w );

[R_global, t_global] = triOptim.optimizeGlobal_Ort_3D( R_c_s_w, t_3D_w );

return
% Code below is deprecated!

%% Final optimization
% Plot map of sampled regions
N = [co.R_c_w];
mask_rot = [co.thereis_line];
N = N(:,mask_rot);
L = cell2mat([co.l]);
% Represent L 2D vectors in sphere normals tangent space
showSamplingSphere(N,L)

%% Set input variables for optimization process
co0 = co;
% precondwithRotDist
precondwithErrFun

%% Optimize rotations
[ R_c_s_w, cov_w, cov_eps_w, err_w, ~, ~ ] = optimRotation( rot_input, R0, true, all_label );
[ R_c_s_nw, cov_nw, cov_eps_nw, err_nw, ~, ~ ] = optimRotation( rot_input, R0, false, all_label );

fprintf('Optimized (weighted) solution R:\n')
disp(R_c_s_w)
fprintf('Optimized (non weighted) solution R:\n')
disp(R_c_s_nw)
fprintf('Angular distance (w and nw): %f\n',angularDistance(R_c_s_w,R_c_s_nw))
if WITHGT
    fprintf('\nGroundtruth value R:\n')
    disp(gt.R_c_s)
end

%% Optimize translations
% TODO: Solve translation fails when triples are not complete
% solveTranslation
% solveTranslation_3D

%% Plot some results
fprintf('Cov of weighted: %f\n', max(eig(cov_w)));
disp( cov_w );
fprintf('Cov of non weighted: %f\n', max(eig(cov_nw)));
disp (cov_nw );

%% Save calibration results
[path_datasets, datasetname]  = fileparts( path );
fbase = fullfile(path_datasets, 'Results');
fext = strcat(typeOfSource,'_',datasetname,...
    stereoLabel,hokuyoLabel,'_',datestr(now,'mm_dd_HH_MM'));
save( strcat(fbase,'R_c_s_W'), 'R_c_s_w', '-ascii' )
save( strcat(fbase,'R_c_s_NW'), 'R_c_s_nw', '-ascii' )
% save( fullfile(path, 't_c_s_w'), 't_c_s_w', '-ascii' )
% save( fullfile(path, 't_c_s_nw'), 't_c_s_nw', '-ascii' )
% saveConfigFile( fullfile(fbase,strcat('W_',fext,'.out')),...
%     struct('R',R_c_s_w,'cov',cov_w,'cov_eps',cov_eps_w ));
% saveConfigFile( fullfile(fbase,strcat('NW_',fext,'.out')),...
%     struct('R',R_c_s_nw,'cov',cov_nw,'cov_eps',cov_eps_nw ));
indexes = [imgs.file_idx];
saveConfigFile( fullfile(fbase,strcat('W_',fext,'.out')),...
    R_c_s_w, cov_w, cov_eps_w, datasetname, stereoLabel, hokuyoLabel, indexes );
saveConfigFile( fullfile(fbase,strcat('NW_',fext,'.out')),...
    R_c_s_nw,cov_nw,cov_eps_nw, datasetname, stereoLabel, hokuyoLabel, indexes );

% plotCalibration( imgs, scans, gt );

% s7_checkResults
