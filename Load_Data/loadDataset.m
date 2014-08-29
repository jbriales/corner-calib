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
    FID = fopen( fullfile(pwd,'cache','path'), 'w' );
    if ( FID >= 0 )
        fprintf( FID, '%s', path );
        fclose( FID );
    else
        error('Could not write file: %s',fullfile(pwd,'cache','path'))
    end
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
    if isempty( [imgs.ts] )
        warning('Wrong image name format: %s\n Check stereo config',...
                imgFormat);
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

out = struct( );
out.imgs = imgs;
out.scans = scans;
out.imgtrack = imgtrack;
out.scantrack = scantrack;
out.R_c_w = R_c_w;
out.signOfAxis = signOfAxis;
out.path = path;
out.hWin = hWin;
out.hLidar = hLidar;
out.hImg = hImg;
out.WITHGT = WITHGT;
out.gt = gt;
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
