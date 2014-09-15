% Main script for the synthetic simulations with variable noise

% Simulated data for the Extrinsic Calibration of a 2D Lidar and a
% Monocular Camera based on Corner Structures without Pattern

% clear classes
clear; clc;

% Main options:
main_sim_file = fullfile( pwd, 'main_Noise.ini' );
mainOpts = readConfigFile( main_sim_file );
extractStructFields( mainOpts );
clear mainOpts

% Set Rig properties
rig_config_file = fullfile( pwd, 'rig.ini' );
rigOpts = readConfigFile( rig_config_file );
extractStructFields( rigOpts );
clear rigOpts
Rig = CSimRig( eye(3), zeros(3,1), R_c_s, t_c_s,... % Extrinsic options
               N, FOVd, scan_sd, d_range,... % Lidar options
               K, res, f, cam_sd ); % Camera options       

% trihedron = CTrihedron( LPattern );
% trihedron = CTrihedron( LPattern, eye(3), 3*[-1 -1 0]' );
trihedron = CTrihedron( LPattern, eye(3), 0*[0 0 -1]' );
corner = CCorner( expmap( [-1 +1 0], deg2rad(-45) ) );
checkerboard = CCheckerboard( RotationZ(deg2rad(45))*RotationY(deg2rad(45)) );
pattern = { trihedron, corner, checkerboard };

% Generate Rig (Camera) poses
% [R_w_c, t_w_c] = generate_random_poses( );
gen_config_file = fullfile( pwd, 'pose_gen.ini' );
[R_w_Cam, R_w_LRF, t_w_Rig, rand_ang_z, rand_ang_x] = generate_random_poses( gen_config_file, Rig ); % For debug purposes only
rand_ang_x = rad2deg( rand_ang_x );
rand_ang_z = rad2deg( rand_ang_z );
Nsamples = length(t_w_Rig);
corresp  = cell(2,3,Nsamples);

% Set Simulation properties
noise_config_file = fullfile( pwd, 'sim_noise.ini' );
noiseOpts = readConfigFile( noise_config_file );
extractStructFields( noiseOpts );
clear noiseOpts
cam_step = (cam_sd_range(2)-cam_sd_range(1))/(cam_sd_N-1);
cam_sd_n = cam_sd_range(1):cam_step:cam_sd_range(2);
scan_step = (scan_sd_range(2)-scan_sd_range(1))/(scan_sd_N-1);
scan_sd_n = scan_sd_range(1):scan_step:scan_sd_range(2);

% Create the output structure
comp = CComparisonNoise( Nsim, cam_sd_N, scan_sd_N, Nsamples );

tic
optim_config_file = fullfile( pwd, 'optim_config.ini' );
optimOpts = readConfigFile( optim_config_file );
extractStructFields( optimOpts );
clear optimOpts
triOptim = CTrihedronOptimization( K,...
    RANSAC_Rotation_threshold,...
    RANSAC_Translation_threshold,...
    debug_level, maxIters,...
    minParamChange, minErrorChange);
cornerOptim = CCornerOptimization( K,...
    debug_level, maxIters,...
    minParamChange, minErrorChange);

% Start the bucle
for j=1:cam_sd_N
for k=1:scan_sd_N
    % Updates the rig (TODO: change the GT for each simulation?¿?¿)
    Rig = CSimRig( eye(3), zeros(3,1), R_c_s, t_c_s,... % Extrinsic options
               N, FOVd, scan_sd_n(k), d_range,... % Lidar options
               K, res, f, cam_sd_n(j) ); % Camera options   
    for m=1:Nsim

        % Store the observations
        for i=1:Nsamples
                % Correspondences for Kwak's algorithm
                if WITHCORNER
                    % Update reference (LRF) pose in Rig for Corner
                    Rig.updateLRFPose( R_w_LRF{i}, t_w_Rig{i} );
                    corr_ = corner.getCorrespondence(Rig);
                    cornerOptim.stackObservation( corr_ );
                end

                % Correspondences for Vasconcelos and Zhang's algorithm
                % Update reference (LRF) pose in Rig for Checkerboard
                if WITHZHANG
                    Rig.updateLRFPose( R_w_LRF{i}, t_w_Rig{i} );
                    check_corresp{1,i} = checkerboard.p2D; 
                    check_corresp{2,i} = checkerboard.getProjection( Rig.Camera );    
                    check_corresp{3,i} = 1000 * cell2mat(checkerboard.getScan( Rig.Lidar ));
                end

                % Correspondences for trihedron
                % Update reference (Camera) pose in Rig for Trihedron
                if WITHTRIHEDRON
                    Rig.updateCamPose( R_w_Cam{i}, t_w_Rig{i} );
                    co_ = trihedron.getCorrespondence( Rig );
                    triOptim.stackObservation( co_ );
                end
        end

    % Trihedron optimization
    if WITHTRIHEDRON
        triOptim.setInitialRotation( [ 0 -1  0
                                       0  0 -1
                                       1  0  0 ] ); % Updated in RANSAC
        comp.TrihedronComp.optim( triOptim, Rig, m, j, k, WITHRANSAC );
    end
    
    % Kwak optimization
    if WITHCORNER
        % Generate random input (near GT)
        R_aux = Rig.R_c_s + randn(3,3)*0.08;
        [U,S,V] = svd(R_aux);
        Rt0 = [ U*V' , Rig.t_c_s + 0.05*randn(3,1) ];
        
        cornerOptim.setInitialRotation( Rt0(1:3,1:3) );
        cornerOptim.setInitialTranslation( Rt0(1:3,4) );
        
        comp.KwakComp.optim( cornerOptim, m, j, k );    
    end
    
    % Zhang and Vasconcelos optimization
    if WITHZHANG
        [T_planes,lidar_points] = checkerboard.getCalibPlanes( Rig, check_corresp );
        comp.VasconcelosComp.optim( T_planes, lidar_points, m, j, k);
        comp.ZhangComp.optim( T_planes, lidar_points, m, j, k);
    end

end
end
end

toc