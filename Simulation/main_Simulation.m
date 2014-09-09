% Main script for the synthetic simulations

% Simulated data for the Extrinsic Calibration of a 2D Lidar and a
% Monocular Camera based on Corner Structures without Pattern

% Previous assertions:
clear; clc;
close all;
dbstop if error;
kbhit('stop');
kbhit('init');      % For stopping anywhere after keyboard hit

% Control variables:
debug = 0;
signOfAxis = [1 1 -1]';
signOfCam  = [1 1 +1]';

% Generation of synthetic data 
fprintf('Generating synthetic data...');
tic
[img, scan, gt] = generate_random_poses_BAK('N',100,'sd_pixel',0,'sd_laser',0,'outlier_k',0,'withPlot',false);
img             = computeCameraCorner( img );
fprintf(' done in %f seconds.\n', toc);

%% FOR LOOP BEGIN
eframe = size(img.pts,3); %eframe = 2;
for nframe = 1:eframe
    %fprintf('nframe = %d\n',nframe)
    %fprintf('-------------\n')
    %% Go from image to Camera-World data
    % Normals to world planes, reprojection plane directions
    %tic;
    [NL, c, A_co, L_P2] = imgParams2calibratedData( img.x(:,nframe), img.K, img.A_x(:,:,nframe) );
    N_repr = L_P2 ./ repmat( sqrt(sum(L_P2.^2,1)), 3,1 );
    %fprintf('CORNER_CALIB TIME: %f\n',toc);
        
    imgs(nframe).NL = NL;
    imgs(nframe).c  = c;
    imgs(nframe).L  = N_repr;
    
    % Solve world plane normals seen from camera and uncertainty
    R0 = gt.R_w_c(:,:,nframe)';
    %tic
    [R_c_w, A_N, A_eps] = getWorldNormals( R0, NL, c, A_co );
    %fprintf('R_c_w optimization TIME: %f\n',toc);
    imgs(nframe).R_c_w    = R_c_w;
    imgs(nframe).A_R_c_w  = A_N;
    %% Go from LIDAR information to Corner Observation
    %tic
    scans(nframe).xy = scan.pts(:,:,nframe);
    for i = 1:3
        if length(scan.inliers{i,nframe}) < 2
            scan.inliers{i,nframe} = [];   
        end
    end
    scan_sigma = 0.03; % TODO: sigma from .ini or generated
    [l,A_l,A_lh,p,A_p,q,A_q, lin,seg] = computeScanCorner( scans(nframe), scan_sigma, scan.inliers(:,nframe), debug );
    thereis_line   = cellfun(@(x)~isempty(x), l);
    thereis_corner = cellfun(@(x)~isempty(x), q);
    if all(thereis_corner) && all(thereis_line)
        CompleteCO = true;
    else
        CompleteCO = false;
    end
    
    scans(nframe).l = l;
    scans(nframe).A_l = A_l;
    scans(nframe).A_lh = A_lh;
    scans(nframe).p = p;
    scans(nframe).A_p = A_p;
    scans(nframe).q = q;
    scans(nframe).A_q = A_q;
    
    % Update calibration closed estimation
    if CompleteCO % Currently deactivated
        % TODO: calculate the sign 
        [R_w_s, t_w_s] = co2LidarPose( cell2mat(q), signOfAxis );
        R_w_s = abs(R_w_s) .* sign(gt.R_w_s(:,:,nframe));
        t_w_s = abs(t_w_s) .* sign(gt.t_w_s(:,nframe));
        scans(nframe).R_w_s = R_w_s;
        scans(nframe).t_w_s = t_w_s;
        
        % Closed estimation of R_c_s according to current frame
        co(nframe).R_c_s = imgs(nframe).R_c_w * scans(nframe).R_w_s;

    end
    %fprintf('SCAN2CO TIME: %f\n',toc);
    %% Store data for final optimisation
    %tic
    lab = [1 2 3];
    co(nframe).complete = CompleteCO;
    co(nframe).thereis_line = thereis_line;        % Line - rotation correspondence
    co(nframe).lab_line = lab(thereis_line);
    co(nframe).thereis_corner = thereis_corner;    % Corner - translation correspondence
    co(nframe).lab_corner = lab(thereis_corner);

    % For rotation optimization
    co(nframe).R_c_w   = imgs(nframe).R_c_w;     % Normal vectors to world planes (in Camera SR)
    co(nframe).A_R_c_w = imgs(nframe).A_R_c_w;
    co(nframe).l       = scans(nframe).l;          % Direction vector of LIDAR segments (in LIDAR SR)
    co(nframe).A_l     = scans(nframe).A_l;
    co(nframe).A_lh    = scans(nframe).A_lh;

    % For translation optimization
    co(nframe).L_P2    = L_P2;
    co(nframe).N_repr  = N_repr;                   % Normal vectors to reprojection planes (in Camera SR)
    co(nframe).q       = q;                        % 2D (XY) Corner points (in LIDAR SR)
    co(nframe).A_q     = A_q;
    %fprintf('STORAGE TIME: %f\n',toc);
        
end

co0 = co;
s_FinalOptimization
% TODO: Doesn't work
solveTranslation
solveTranslation_3D

