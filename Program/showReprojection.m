% All scripts

% Previous assertions:
% clear corresp corresp_t
clear
dbstop if error
kbhit('stop')
kbhit('init') % For stopping anywhere after keyboard hit

% Control variables:
global WITH_MONTECARLO
WITH_MONTECARLO = false;

SIMULATION = false;
RAWLOG = false;
typeOfSource = 'Rawlog';
typeOfLineDetection = 'Manual';
USER_INITIALISE = true;
if ~exist('USER_INITIALISE','var')
    USER_INITIALISE = false;
end
dataset_Nr = 1;


% Automated dataset selection
switch typeOfSource
    case 'Blender'
        path = fullfile(dropbox,'Blender files/Corner/files');
%         path = fullfile(dropbox,'Blender files/Corner_3/files');
        decimation = 2;
        imgFormat = '%f.png';
        THRES_PRE = 0.2;
        WITHGT = true;
        debugControl = false;
        
        signOfAxis = [1 1 -1]';
        signOfCam  = [1 1 +1]';
        
        LIDAR_tag = [];
        
    case 'Rawlog'
        decimation = 1;
        img_decimation = 25;
        THRES_PRE = 0.02;
        WITHGT = false;
        debugControl = false;
        
        signOfAxis = [1 1 +1]';
        signOfCam  = [1 1 -1]';

%         C_imgFormat = {'_%f.jpg','_%f.jpg','CAM_%f.jpg','CAM_%f.jpg','CAM_%f.jpg','CAM_%f.jpg'};
        str = strtrim( evalc('system(''hostname'');') );
        switch str
            case 'jesus-desktop'
                 C_LIDAR_tag = {'_LASER_HOKUYO_UTM', '_LASER_HOKUYO_UTM', '_LASER_HOKUYO', '_LASER_HOKUYO_A', '_LASER_HOKUYO_A', '_LASER_HOKUYO'};
%                  path = '/media/storage/Research/Corner_Datasets/Corner';
                 path = '/media/storage/Research/Corner_Datasets/Check';
                 C_path = { '1', '2', '3', '4', '5', '6' };
            case 'jesus-P5QL-PRO'
                C_LIDAR_tag = {'_LASER_HOKUYO_A'};
                path = '/media/cloud/Datasets/Corner';
                C_path = { 'A' };
            otherwise
                error('Choose automated dataset location')
        end
%         imgFormat = C_imgFormat{dataset_Nr};
        imgFormat = 'CAM_%f.jpg';
        LIDAR_tag = C_LIDAR_tag{dataset_Nr};
        path = strcat( path, C_path{dataset_Nr} );
        clear str
end

if ~exist('path','var')
    path = uigetdir('','Choose dataset parent folder');
end

if ~exist( fullfile( path, 'user_complete_initialization.mat' ), 'file' )
    USER_INITIALISE = true;
end
if ~USER_INITIALISE
    tic
    load('user_complete_initialization');
    toc
    
% %     hDocked = figure('Name',['Calibration in ',path], 'WindowStyle','docked');
%     hDocked = figure('Name',['Calibration in ',path]);
%     hImageTrack = subplot(1,2,2);
%     imshow( rgb2gray(imgs(1).I) )
% 
%     hLidarPlane = subplot(121);
%     hold on, title('First scan')
% %     [lin, seg] = manualSetLidarSegments( scans(1).xy, lin, seg );
%     [lin] = manualSetLidarSegments( scans(1).xy, lin );
    
else
    %% Set initial data and user input
    blur = false;
    imgs  = loadCamData( path, imgFormat, imgFormat, blur, img_decimation );
    % imgsStack = stackImages( imgs );
    % implay( imgsStack, 80 )
    % pause
    scans = readLidarData(typeOfSource, path, LIDAR_tag);
    
    % Get indexes for timestamp limits
    time = [scans.ts];
    im_tspan = [imgs.ts];
    scan_first_idx = find( time>=im_tspan(1), 1 );
    scan_last_idx  = find( time>=im_tspan(end), 1 );
    
    scans(scan_last_idx+1:end) = [];
    scans(1:scan_first_idx-1) = [];
    clear time im_tspan scan_first_idx scan_last_idx
    
    % vnframe = 1:length(imgs);
    % decimation = 5;
    if decimation > 1
        scans = scans(1:decimation:end);
    end
    vnframe = 1:length(scans);
    
    save(fullfile( path, 'user_complete_initialization' ))
end


%% Show reprojection
%% Begin for loop for data extraction
% co = repmat( struct('cam',[], 'scan',[]), length(imgs), 1 );
CompleteCO = true;
repr_planes = [];
scan_points = [];
im_frame = 0;
CHECK_IMAGE = false;
nco = 0;

% Plot reprojection
figure('Position',[687   366   560   420]), hold on
% plotframe( eye(4), 0.1, 'C', 'k' )
plotframe( eye(4), 0.1, 'S', 'k' )

load( fullfile(path, 'R_c_s_w') )
load( fullfile(path, 't_c_s_w') )
R_c_s = R_c_s_w;
% t_c_s = t_c_s_w;
t_c_s = [ 0.005 -0.075 -0.070 ]';
Rt_c_s = [R_c_s t_c_s];
T_c_s = [Rt_c_s ; 0 0 0 1];

% plotframe( T_c_s, 0.1, 'S', 'b' )
plotframe( inv(T_c_s), 0.1, 'C', 'b' )
axis equal, rotate3d on
keyboard
close

hFig = figure('Position',[687   366   560   420]);
title('Image with superimposed x-coloured LIDAR points')
for nframe=vnframe
    % Set algorithm input
    scan = scans(nframe);
    
    if im_frame < find( [imgs.ts] <= scan.ts, 1, 'last' );
        new_image_reached = true;
    else
        new_image_reached = false;
    end
    
    if new_image_reached
        nco = nco + 1; % If new image is reached new Corner Observation will be added
%         img_occ = cellfun(@(x)~isempty(x),{imgs.ts}, 'UniformOutput', false);
%         img_occ = [img_occ{:}];

        im_frame = find( [imgs.ts] <= scan.ts, 1, 'last' );
        img = imgs(im_frame);

        cla
        reprojectScan2Img( img.I, scan.xy, img.K, Rt_c_s, hFig )
%         pause(0.05)
        keyboard
    end
    
end
warning('Reprojection plot ended!')