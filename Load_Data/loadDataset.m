function out = loadDataset( userOpts )
% out = loadDataset( userOpts )
% Receives the fields specified in userOpts and load data related to
% dataset accordingly
% The output data is stored in a struct
% 
% See also extractStructFields
      
% Copy struct values to variables with field names
extractStructFields( userOpts );

if exist(fullfile(pwd,'cache','path'),'file')
    % Read path value from cache file
    FID = fopen( fullfile(pwd,'cache','path'), 'r' );
    if ( FID >= 0 )
        path = fgetl( FID );
        fclose(FID);
    else
        error('Could not read file: %s',fullfile(pwd,'cache','path'))
    end
    % Look chosen dataset in parent folder
    hint_path = fullfile( fileparts(path), datasetname );
    if ~exist(hint_path,'dir')
        % If do not exist: uigetdir in parent folder
        path = uigetdir(fileparts(path),'Choose dataset parent folder');
    else
        path = hint_path;
    end
else
    % If there is no path stored in cache: uigetdir
    path = uigetdir(pwd,'Choose dataset parent folder');
end

if USER_INITIALISE
    if ~exist(fullfile(pwd,'cache'),'dir')
        mkdir(fullfile(pwd,'cache'));
        warning('cache dir does not exist in %s.\nIt has been created.',pwd);
    end
    FID = fopen( fullfile(pwd,'cache','path'), 'w' );
    if ( FID >= 0 )
        fprintf( FID, '%s', path );
        fclose( FID );
    else
        error('Could not write file: %s',fullfile(pwd,'cache','path'))
    end
end

% Assure auxiliar folders exist in dataset folder
aux_folders = {'cache','meta_img','meta_laser'};
for fold = aux_folders
    fold = fold{:};
    if ~exist(fullfile(path,fold),'dir')
        mkdir(path,fold)
    end
end

% Set plotting space
if exist(fullfile(pwd,'cache','winPosition.mat'),'file')
    load( fullfile(pwd,'cache','winPosition.mat'), 'winPosition' );
    hWin = figure('Name',['Calibration in ',path],...
                     'WindowStyle',WindowStyle,...
                     'Position', winPosition,...
                     'Visible', win_visibility);
else
    hWin = figure('Name',['Calibration in ',path],...
                     'WindowStyle',WindowStyle,...
                     'Visible', win_visibility);
end
if videoFigure
%     hLidar = subplot(2,2,4);
%     hLidar = subplot_tight(3,3,[3 6]);
    hLidar = subplot(6,3,[3 6 9]);
    title('Camera capture')
%     hImg   = subplot(2,2,[1]);
%     hImg   = subplot_tight(3,3,[1 2 4 5]);
    hImg   = subplot(6,3,[1 2 4 5 7 8 10 11]);
    
    hVid   = subplot(6,3,[12 15 18]);
    title('LRF 2D scan')
    set(gcf,'color','w');
else
    hLidar = subplot(1,2,1);
    hImg = subplot(1,2,2);
end
% Roll LRF axis
% subplot( hLidar );
% camroll(90);

if ~USER_INITIALISE
    tic
    load( fullfile(path,'cache','initialisation') );
    toc
    
    % hDocked = figure('Name',['Calibration in ',path], 'WindowStyle','docked');
    % hDocked = figure('Name',['Calibration in ',path], 'Position', [350 345 1108 571]);
    subplot( hImg );
    imshow( rgb2gray(imgs(1).I) )
    
    subplot( hLidar );
    hold on, title('Current scan segmentation')
    plot( scans(1).xy(1,:), scans(1).xy(2,:), '.k' );
    axis equal
    plotLIDARframe();
else
    %% AUTOMATED DATASET SELECTION AND CONFIG
    % TODO: Change signOfAxis and signOfCam by automated
    switch typeOfSource
        case 'Blender'
            imgFormat = 'CAM_%f.png';
            WITHGT = true;
            gt = loadCalibrationGT();
            signOfAxis = [1 1 -1]';
            signOfCam  = [1 1 +1]';
            LIDAR_tag = [];
        case 'Rawlog'
            WITHGT = false;
            gt = [];
            signOfAxis = [1 1 +1]';
            signOfCam  = [1 1 -1]';
            %         str = strtrim( evalc('system(''hostname'');') );
            imgFormat = strcat('CAM',stereoLabel,'_%f.jpg');
            LIDAR_tag = strcat(LIDAR_tag,hokuyoLabel);
    end

    %% LOAD INITIAL DATA AND USER INPUT
    SConfigCam = readConfigFile( fullfile(path,'rig.ini'),'[Cam]' );
    Cam = CRealCamera( SConfigCam );
    SConfigLRF = readConfigFile( fullfile(path,'rig.ini'),'[LRF]' );
    LRF = CRealLidar( SConfigLRF );
    % Create simulated object with given configuration:
    CConfigCam = struct2cell( SConfigCam );
    CConfigLRF = struct2cell( SConfigLRF );
    SimRig = CSimRig( eye(3), zeros(3,1),... % Extrinsic options
               eye(3), zeros(3,1),...
               CConfigLRF{:},... % Lidar options
               CConfigCam{:} ); % Camera options 
    
    % Set handle for figure
    Cam.hFig = hImg;
    LRF.hFig = hLidar;
    subplot( LRF.hFig );
    camroll(90);
    
    % Load images
    blur = false;
    imgs = loadCamData( path, imgFormat, blur, decimation, img_idxs );
    if isempty( imgs )
        warning('No images selected. The program will finish')
        return
    end
    if isempty( [imgs.ts] )
        warning('Wrong image name format: %s\n Check stereo config',...
                imgFormat);
        return
    end
    
%     Temporarily unavailable
%     if VIDEO_PREVIEW
%         imgsStack = stackImages( imgs );
%         implay( imgsStack, 80 )
%         pause
%         clear imgsStack
%     end
    
    % Load LIDAR
%     scans = loadLidarData(typeOfSource, path, LIDAR_tag);
    rawlogLRF = CRawlogLRF( fullfile(path,'laser') );
    scans = rawlogLRF.scans;
    
    % TODO: Improve LIDAR-cam synchronisation
    % Remove images before and after scan limits:
%     plotTs( scans, imgs )
    imgs  = cropImgsArray( scans, imgs );
    
    % Plot scans and images timestamps (prior):
    if WITH_plotTs
        figure('Name','timestamp LIDAR and Camera'), hold on
        plotTs( scans, imgs, '.k', 'ok' );
    end
    % Remove scans before and after images limits:
%     plotTs( scans, imgs )
    %scans = cropScansArray( scans, imgs );
    % Set image-scan synchronised pairs by finding scan closest to each image and assign:
%     plotTs( scans, imgs )
    [scans, imgs] = filterClosest( scans, imgs, syncr_gap );
%     [imgs, delta_ts] = filterClosestFrames( scans, imgs );
%     [scans, delta_ts] = filterClosestScans( scans, imgs );
    
    % Plot scans and images timestamps (after correspondences):
    if WITH_plotTs
        plotTs( scans, imgs, '+r','sr' );
        pause
        close
    end
    
%     scans  = CScan( scans ); % Not needed if using CRawlogLRF
    frames = CFrame( imgs );
% %     % User input: Set initial points in image for tracking
% % %     Cam.setFrame( frames(1) );
% %     if exist(imgs(1).metafile,'file')
% % %         load( imgs(1).metafile, '-mat', 'imgtrack' );
% %         load( fullfile(imgs(1).path,'meta_img',imgs(1).metafile),...
% %             '-mat', 'imgtrack' );
% %     else
% % %         imgs(1).I = imread( imgs(1).path );
% %         imgs(1).I = imread( fullfile(imgs(1).path,'img',imgs(1).file) );
% %         imgtrack = initialisation_calib( rgb2gray(imgs(1).I) );
% %     end
% %     
% % %     imgs(1).I = imread( imgs(1).path );
% %     imgs(1).I = imread( fullfile(imgs(1).path,'img',imgs(1).file) );
% %     subplot(hImg);
% %     imshow( rgb2gray(imgs(1).I) );
    
    % User input: Set segments of interest in scan observation
% %     scan = scans(1);
% %     switch typeOfSource
% %         case 'Blender'
% %             scan.xy = loadBlenderPCD( scan.path );
% %         case 'Rawlog'
% %     end
% %     hLidar = subplot(121);
% %     hold on, title('First scan')
% %     if exist(scan.metafile,'file')
% %         load( scan.metafile, '-mat', 'scantrack' );
% %         plot( scan.xy(1,:), scan.xy(2,:), '.k' )
% %         axis equal
% %         plotLIDARframe
% %     else
% % %         scantrack = manualSetScanlines( scans(1).xy, [1 1 1] );
% %         scantrack = manualSetScanlines( scan.xy, [1 1 1], scan.metafile );
% %     end
    
%     [R_w_c, t_w_c] = manualSetCameraPose( signOfCam, 0 );
    % [R_w_c, t_w_c] = manualSetCameraPose( );
        
    % Set rotation R_c_w initial estimate
%     R_c_w = R_w_c'; % Inverse of R_w_c
    
%     save( fullfile(path,'cache','initialisation'),...
%           'imgs', 'scans', 'imgtrack', 'scantrack', 'R_c_w', 'WITHGT',...
%           'signOfAxis', 'gt')
end
winPosition = get(gcf,'Position');
save( fullfile(pwd,'cache','winPosition'), 'winPosition' )

out = struct( );
out.imgs = imgs;
out.scans = scans;
% out.imgtrack = imgtrack;
% out.scantrack = scantrack;
% out.R_c_w = R_c_w;
% out.signOfAxis = signOfAxis;
out.path = path;
out.hWin = hWin;
out.hLidar = hLidar;
out.hImg = hImg;
out.hVid = hVid;
out.WITHGT = WITHGT;
out.gt = gt;
out.Cam = Cam;
out.LRF = LRF;
out.SimRig = SimRig;
out.frames = frames;
end

%% Local functions
function imgs = cropImgsArray( scans, imgs )
scan_ts = [scans.ts];
cam_ts  = [imgs.ts];
cam_first_idx = find( cam_ts>=scan_ts(1), 1 );
cam_last_idx  = find( cam_ts>=scan_ts(end), 1 );

imgs(cam_last_idx+1:end) = [];
imgs(1:cam_first_idx-1) = [];
end

function scans = cropScansArray( scans, imgs )
scan_ts = [scans.ts];
cam_ts  = [imgs.ts];
scan_first_idx = find( scan_ts>=cam_ts(1), 1 );
scan_last_idx  = find( scan_ts>=cam_ts(end), 1 );

scans(scan_last_idx+1:end) = [];
scans(1:scan_first_idx-1) = [];
end

function [scans_, imgs_] = filterClosest( scans, imgs, threshold )
scan_ts = [scans.ts];
counter = 1;
% threshold = 1;
for i=1:length(imgs)
    diff = scan_ts - imgs(i).ts;
    [delta_ts,I] = min( abs(diff) );
    
    if delta_ts < threshold
        scans_(counter) = scans(I);
        scans_(counter).delta_ts = delta_ts;
        imgs_(counter)  = imgs(i);
        counter = counter + 1;
    end
end
end

function [scans, delta_ts] = filterClosestScans( scans, imgs )
scan_ts = [scans.ts];
delta_ts = zeros(1,length(imgs));
for i=1:length(imgs)
    diff = scan_ts - imgs(i).ts;
    [delta_ts(i),I] = min( abs(diff) );
    scans_(i) = scans(I);
end
for i=1:length(imgs)
    scans_(i).delta_ts = delta_ts(i);
end
scans = scans_; % Substitute array
clear scans_
end

function [imgs, delta_ts] = filterClosestFrames( scans, imgs )
imgs_ts = [imgs.ts];
% scan_ts = [scans.ts];
delta_ts = zeros(1,length(imgs));
for i=1:length(scans)
    diff = scans(i).ts - imgs_ts;
    [delta_ts(i),I] = min( abs(diff) );
    imgs_(i) = imgs(I);
end
for i=1:length(scans)
    scans(i).delta_ts = delta_ts(i);
end
imgs = imgs_; % Substitute array
clear imgs_
end

function plotTs( scans, imgs, format_scan, format_img )
scan_ts = [scans.ts];
cam_ts = [imgs.ts];
time_min = min( [scan_ts cam_ts] );

plot( scan_ts-time_min, ones(1,length( scan_ts )), format_scan )
plot( cam_ts-time_min, ones(1,length( cam_ts )), format_img )
legend('LIDAR','Camera');
if isfield(scans,'delta_ts')
    title(sprintf('Max ts dist: %f',max([scans.delta_ts])));
end
% pause
% close
end
