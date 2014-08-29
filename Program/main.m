% All scripts

% Previous assertions:
clear
close all
dbstop if error % Set debug mode on if error occurs
kbhit('stop')
kbhit('init') % Track keyboard hits

set(0,'defaultfigureposition', [642   503   560   420])

% Control variables:
global WITH_MONTECARLO
WITH_MONTECARLO = false;

userOpts = readConfigFile( fullfile(pwd,'main.ini') );
extractStructFields( userOpts );

all_dataset_data = loadDataset( userOpts );
extractStructFields( all_dataset_data );

debug = 0;
[img_params, A_img_params, imgtrack, checkImage] = ...
    corner_calib(imgtrack, rgb2gray(imgs(1).I), debug);

CompleteCO = true;
repr_planes = [];
scan_points = [];
im_frame = 0;
CHECK_IMAGE = false;
nco = 0;

%% FOR LOOP BEGIN
for nframe=1:length(scans)
    fprintf('nframe = %d\n',nframe)
    fprintf('-------------\n')
    figure( hWin )
    
    if upper(kbhit) == 'D' % To stop after pressing 'L' key
        disp('D was hit: Debug activation')
        keyboard
    end
    
    % Set algorithm input
    scan = scans(nframe);
    
    if im_frame < find( [imgs.ts] <= scan.ts, 1, 'last' );
        new_image_reached = true;
    else
        new_image_reached = false;
    end
    
    if new_image_reached
        nco = nco + 1; % If new image is reached new Corner Observation will be added
        
        im_frame = find( [imgs.ts] <= scan.ts, 1, 'last' );
        img = imgs(im_frame);
        
        if WITHGT
            gt = loadGT( gt, path, im_frame );
            co(nco).gt = gt;
        end
        
        %% Go from image to Camera-World data
        imgtrack0 = imgtrack;
        img_gray = rgb2gray( img.I );
        
        subplot( hImg )
        fprintf('Image tracking in frame %d\n',nframe)
        tic
        debug = 0;
        [img_params, A_img_params, imgtrack, CHECK_IMAGE] = ...
            corner_calib(imgtrack, img_gray, debug);
        fprintf('CORNER_CALIB TIME: %f\n',toc)
        clear debug
        
        if CHECK_IMAGE
            subplot( hImg )
            corner_calib_plot(imgtrack.x, imgtrack0.x, img_gray)
            
            imgtrack = initialisation_calib( img_gray );

            [~, ~, imgtrack, ~] = corner_calib(imgtrack, img_gray, debug);
            [img_params, A_img_params, imgtrack, ~] = ...
                corner_calib(imgtrack, img_gray, debug);
            % Twice to avoid gradient direction alarm from initialisation
            CHECK_IMAGE = false;
        end
        x = img_params;
        A_x = A_img_params;
        [NL, c, A_co, L_P2] = imgParams2calibratedData( x, img.K, A_x );
        N_repr = L_P2 ./ repmat( sqrt(sum(L_P2.^2,1)), 3,1 );
        
        %% Plot image tracking
        tic
        subplot( hImg )
        corner_calib_plot(img_params, imgtrack0.x, img_gray)
        fprintf('IMAGE PLOT TIME: %f\n',toc)
            
        % Solve world plane normals seen from camera and uncertainty
        R0 = R_c_w;
        tic
        [R_c_w, A_R_c_w, A_eps_c_w] = getWorldNormals( R0, NL, c, A_co );
        fprintf('R_c_w optimization TIME: %f\n',toc)
        if WITHGT
            disp('Compare computed R_c_w to GT')
            disp(R_c_w), disp(gt.R_c_w)
            fprintf('Angular distance with GT rotation: %f\n', angularDistance(R_c_w, gt.R_c_w))
            if angularDistance(R_c_w, gt.R_c_w) > 3
                CHECK_IMAGE = true;
            end
        end
        
        % Check displacement from previous frame
        if exist('R_c_w0','var')
            fprintf('Angular distance from image %d to %d: %e\n', im_frame-1, im_frame, angularDistance(R_c_w0,R_c_w))
        end
        R_c_w0 = R_c_w; % Store previous for comparison       
                
        debug_im_reprojection = false;
        if debug_im_reprojection
            % Check world projection with obtained rotation
            figure
            imshow( img.I ); hold on
            checkRotationReprojection( c, R_c_w, img.K )
            pause, close
        end
    end

    %% Go from LIDAR information to Corner Observation
    % User reallocation of LIDAR segments
    if upper(kbhit) == 'L' % To stop after pressing 'L' key
        disp('L was hit: Manual reallocation of segments')
        keyboard
        
        if ~exist('thereis_line','var') || all(thereis_line)
            selection_mask = input('Input segment selection mask [x y z]: ');
        else
            selection_mask = ~thereis_line;
        end
        
        subplot( hLidar )
        scantrack_aux = manualSetScanlines( scan.xy, selection_mask );
        for k=1:3
            if selection_mask(k) % Substitutes only user selected segments
                scantrack(k) = scantrack_aux(k);
            end
        end
    end
        
    debug = 0;
    tic
    [scantrack, inPts, lost] = updateSegments( scan, scantrack, debug );
    [l,A_l,A_lh,p,A_p,q,A_q, lin,seg] = computeScanCorner( scan, {scantrack.all_inliers}, debug );
    [scantrack.lin] = deal( lin{:} );
    [scantrack.seg] = deal( seg{:} );

    fprintf('SCAN2CO TIME: %f\n',toc)
    thereis_line   = cellfun(@(x)~isempty(x), l);
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
    tic
    subplot(hLidar)
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
       
    % Update calibration closed estimation
    tic
    if new_image_reached && CompleteCO % Currently deactivated
        % TODO: signOfAxis could change and this has to be automated or
        % problems can arise (complex numbers, etc.)
        [R_w_s, t_w_s] = co2LidarPose( cell2mat(q), signOfAxis );
        
        % Closed estimation of R_c_s according to current frame
        co(nco).R_c_s = R_c_w * R_w_s;
        
        if 0
            % Average of estimated rotations
            R_c_s = sum( reshape([co.R_c_s],3,3,[]), 3 );
            [U,~,V] = svd( R_c_s );
            R_c_s = U*V';
            
            % Solve t_c_s (closed solution with Lidar and Cam poses)
            % Add new complete observation
            repr_planes = [repr_planes, imgs(im_frame).L]; %#ok<AGROW>
            scan_points = [scan_points, cell2mat(scans(nframe).q)]; %#ok<AGROW>
            b = - dot( repr_planes, R_c_s(:,1:2) * scan_points, 1 )';
            A = repr_planes';
            
            t_c_s = A \ b;
        end
    end
    fprintf('Current frame R_c_s estimation TIME: %f\n',toc)
    
    %% Store data for final optimization
    tic
    if new_image_reached
        lab = [1 2 3];
        co(nco).complete = CompleteCO;
        co(nco).thereis_line = thereis_line; % Line - rotation correspondence
        co(nco).lab_line = lab(thereis_line);
        co(nco).thereis_corner = thereis_corner; % Corner - translation correspondence
        co(nco).lab_corner = lab(thereis_corner);
        
        % For rotation optimization
        co(nco).R_c_w   = R_c_w; % Normal vectors to world planes (in Camera SR)
        co(nco).A_R_c_w = A_R_c_w;
        co(nco).l       = l; % Direction vector of LIDAR segments (in LIDAR SR)
        co(nco).A_l     = A_l;
        
        % For translation optimization
        co(nco).N_repr  = N_repr; % Normal vectors to reprojection planes (in Camera SR)
        co(nco).L_P2    = L_P2;
        co(nco).A_lh    = A_lh;
        co(nco).q       = q; % 2D (XY) Corner points (in LIDAR SR)
        co(nco).A_q     = A_q;
        
        % For final checking and debug
        check(nco).cam_frame  = im_frame;
        check(nco).scan_frame = nframe;
        check(nco).scan_inPts = inPts;
        check(nco).scan_track = scantrack;
        check(nco).img_track  = imgtrack;
    end
    fprintf('STORAGE TIME: %f\n',toc)

end
save( fullfile( path, 'cache', strcat('preoptimization',stereoLabel,hokuyoLabel) ) )
keyboard

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
[ R_c_s_w, err_w, ~, ~ ] = optimRotation( rot_input, R0, true, all_label );
[ R_c_s_nw, err_nw, ~, ~ ] = optimRotation( rot_input, R0, false, all_label );

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
solveTranslation
solveTranslation_3D

%% Save calibration results
[path_datasets, datasetname]  = fileparts( path );
fbase = fullfile(path_datasets, 'Results',...
	strcat(typeOfSource,'_',datasetname,...
    stereoLabel,hokuyoLabel,'_',datestr(now,'mm_dd_HH_MM'),'_'));
save( strcat(fbase,'R_c_s_W'), 'R_c_s_w', '-ascii' )
save( strcat(fbase,'R_c_s_NW'), 'R_c_s_nw', '-ascii' )
% save( fullfile(path, 't_c_s_w'), 't_c_s_w', '-ascii' )
% save( fullfile(path, 't_c_s_nw'), 't_c_s_nw', '-ascii' )

% plotCalibration( imgs, scans, gt );

% s7_checkResults
