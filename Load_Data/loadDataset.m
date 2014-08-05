function [imgs, scans, imgtrack, scantrack,R_c_w,signOfAxis,...
          path, hWin, hLidar, hImg, WITHGT, gt ] = loadDataset( userOpts )

% Copy struct values to variables with field names
vars = fieldnames( userOpts );
for k=1:length(vars)
    field = vars{k};
    eval([ field,'= userOpts.',field,';' ]);
end
clear k field vars userOpts

if exist(fullfile(pwd,'cache','path.mat'),'file')
    load( fullfile(pwd,'cache','path'),'-mat', 'path' );
else
    path = pwd;
end
if USER_INITIALISE
    path = uigetdir(path,'Choose dataset parent folder');
    save( fullfile(pwd,'cache','path'), 'path' );
end

% Set plotting space
if exist(fullfile(pwd,'cache','winPosition.mat'),'file')
    load( fullfile(pwd,'cache','winPosition.mat'), 'winPosition' );
    hWin = figure('Name',['Calibration in ',path],...
                     'WindowStyle',WindowStyle,...
                     'Position', winPosition );
else
    hWin = figure('Name',['Calibration in ',path],...
                     'WindowStyle',WindowStyle);
end
hLidar = subplot(1,2,1);
hImg = subplot(1,2,2);

if ~USER_INITIALISE
    tic
    load( fullfile(path,'cache','initialisation') );
    toc
    
    % hDocked = figure('Name',['Calibration in ',path], 'WindowStyle','docked');
%     hDocked = figure('Name',['Calibration in ',path], 'Position', [350 345 1108 571]);
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
            imgFormat = '%f.png';
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
            LIDAR_tag = strcat('_LASER_HOKUYO',hokuyoLabel);
    end

    %% LOAD INITIAL DATA AND USER INPUT
    % Load images
    blur = false;
    imgs = loadCamData( path, imgFormat, imgFormat, blur, decimation );
    if isempty( imgs )
        warning('No images selected. The program will finish')
        return
    end
    
    if VIDEO_PREVIEW
        imgsStack = stackImages( imgs );
        implay( imgsStack, 80 )
        pause
        clear imgsStack
    end
    
    % Load LIDAR
    scans = loadLidarData(typeOfSource, path, LIDAR_tag);
    
    % TODO: Improve LIDAR-cam synchronisation
    % Remove images before and after scan limits:
%     plotTs( scans, imgs )
    imgs  = cropImgsArray( scans, imgs );
    % Remove scans before and after images limits:
%     plotTs( scans, imgs )
    scans = cropScansArray( scans, imgs );
    % Find scan closest to each image and assign:
%     plotTs( scans, imgs )
    scans = filterClosestScans( scans, imgs );
    % Plot scans and images timestamps:
%     plotTs( scans, imgs )
    
    % User input: Set initial points in image for tracking
    imgtrack = initialisation_calib( rgb2gray(imgs(1).I) );
    
    subplot(hImg);
    imshow( rgb2gray(imgs(1).I) );
    
    % User input: Set segments of interest in scan observation
    hLidar = subplot(121);
    hold on, title('First scan')
    scantrack = manualSetScanlines( scans(1).xy, [1 1 1] );
    
    [R_w_c, t_w_c] = manualSetCameraPose( signOfCam, 0 );
    % [R_w_c, t_w_c] = manualSetCameraPose( );
        
    % Set rotation R_c_w initial estimate
    R_c_w = R_w_c'; % Inverse of R_w_c
    
    save( fullfile(path,'cache','initialisation'),...
          'imgs', 'scans', 'imgtrack', 'scantrack', 'R_c_w', 'WITHGT',...
          'signOfAxis', 'gt')
end
winPosition = get(gcf,'Position');
save( fullfile(pwd,'cache','winPosition'), 'winPosition' )
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

function scans = filterClosestScans( scans, imgs )
scan_ts = [scans.ts];
for i=1:length(imgs)
    diff = scan_ts - imgs(i).ts;
    [~,I] = min( abs(diff) );
    scans_(i) = scans(I);
end
scans = scans_; % Substitute array
clear scans_
end

function plotTs( scans, imgs )
scan_ts = [scans.ts];
cam_ts = [imgs.ts];
time_min = min( [scan_ts cam_ts] );
figure('Name','timestamp LIDAR and Camera'), hold on
plot( scan_ts-time_min, 1, '.b' )
plot( cam_ts-time_min, 1, 'or' )
legend('LIDAR','Camera')
pause
close
end
