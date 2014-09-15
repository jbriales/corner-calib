% Main script for the synthetic simulations with variable noise

% Simulated data for the Extrinsic Calibration of a 2D Lidar and a
% Monocular Camera based on Corner Structures without Pattern

% clear classes
clear; clc;

% Generate Rig (Camera) poses
% [R_w_c, t_w_c] = generate_random_poses( );
gen_config_file = fullfile( pwd, 'pose_gen.ini' );
[R_w_c, t_w_c, rand_ang_z, rand_ang_x] = generate_random_poses( gen_config_file ); % For debug purposes only
rand_ang_x = rad2deg( rand_ang_x );
rand_ang_z = rad2deg( rand_ang_z );
Nsamples = length(R_w_c);
corresp  = cell(2,3,Nsamples);

% Set Rig properties
rig_config_file = fullfile( pwd, 'rig.ini' );
rigOpts = readConfigFile( rig_config_file );
extractStructFields( rigOpts );
clear rigOpts
Rig = CSimRig( eye(3), zeros(3,1), R_c_s, t_c_s,... % Extrinsic options
               N, FOVd, scan_sd, d_range,... % Lidar options
               K, res, f, cam_sd ); % Camera options                

trihedron = CTrihedron( LPattern );
corner = CCorner( expmap( [-1 +1 0], deg2rad(-45) ) );
checkerboard = CCheckerboard( RotationZ(deg2rad(45))*RotationY(deg2rad(45)) );
pattern = { trihedron, corner, checkerboard };

corner_corresp = cell(2,3,Nsamples);

tic
optim_config_file = fullfile( pwd, 'optim_config.ini' );
optimOpts = readConfigFile( optim_config_file );
extractStructFields( optimOpts );
clear optimOpts
triOptim = CTrihedronOptimization( K,...
    RANSAC_Rotation_threshold,...
    RANSAC_Translation_threshold,...
    debug_level, maxIters );

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
        % Update reference (Camera) pose in Rig
        Rig.updatePose( R_w_c{i}, t_w_c{i} );

        % Correspondences for Kwak's algorithm
        corr_ = corner.getCorrespondence(Rig); 
        for p = 1:2
            for q = 1:3
                corner_corresp{p,q,i} = corr_{p,q};
            end
        end   

        % Correspondences for Vasconcelos and Zhang's algorithm
        check_corresp{1,i} = checkerboard.p2D; 
        check_corresp{2,i} = checkerboard.getProjection( Rig.Camera );    
        check_corresp{3,i} = 1000 * cell2mat(checkerboard.getScan( Rig.Lidar ));

        % Correspondences for trihedron 
    %     trihedron.plotScene(Rig.Camera, Rig.Lidar);
    %     close;
        co_ = trihedron.getCorrespondence( Rig );
    %     co(i) = co_;
        triOptim.stackObservation( co_ );
    end

% Trihedron optimization
triOptim.setInitialRotation( [ 0 -1  0
                               0  0 -1
                               1  0  0 ] ); % Updated in RANSAC
comp.TrihedronComp.optim( triOptim, m, j, k );

end
end
end

toc