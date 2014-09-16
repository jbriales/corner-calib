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

% Create patterns
trihedron = CTrihedron( LTrihedron, eye(3), 0*[0 0 -1]' );
corner = CCorner( LCorner, expmap( [-1 +1 0], deg2rad(-45) ) );
checkerboard = CCheckerboard( LCheckerboard, RotationZ(deg2rad(45))*RotationY(deg2rad(45)) );
pattern = { trihedron, corner, checkerboard };

% Set Simulation properties
noise_config_file = fullfile( pwd, 'sim_noise.ini' );
noiseOpts = readConfigFile( noise_config_file );
extractStructFields( noiseOpts );
clear noiseOpts

cam_sd_n  =   logspace(cam_sd_range(1),cam_sd_range(2),cam_sd_N);
scan_sd_n = 3*logspace(scan_sd_range(1),scan_sd_range(2),scan_sd_N);
N_co_N    = size(N_co_n,2);

% Create the output structure for comparison
comp = CComparisonNoise( Nsim, cam_sd_n, scan_sd_n, N_co_n );

% Create objects for storage and optimization
optim_config_file = fullfile( pwd, 'optim_config.ini' );
optimOpts = readConfigFile( optim_config_file );
extractStructFields( optimOpts );
clear optimOpts

% Start the loop
Nobs = N_co_n(end);
for idx_cam_sd=1:cam_sd_N
for idx_scan_sd=1:scan_sd_N
    % Updates the rig
    Rig = CSimRig( eye(3), zeros(3,1), R_c_s, t_c_s,... % Extrinsic options
        N, FOVd, scan_sd_n(idx_scan_sd), d_range,... % Lidar options
        K, res, f, cam_sd_n(idx_cam_sd) ); % Camera options
    fprintf('Simulations for cam_sd = %.2E and scan_sd = %.2E\n',cam_sd_n(idx_cam_sd), scan_sd_n(idx_scan_sd));
    for idx_Nsim=1:Nsim
        % Create optimization objects and random poses
        if WITHTRIHEDRON
            clearvars triOptim;
            triOptim = CTrihedronOptimization( K,...
                RANSAC_Rotation_threshold,...
                RANSAC_Translation_threshold,...
                debug_level, maxIters,...
                minParamChange, minErrorChange);
            
            gen_config_file = fullfile( pwd, 'pose_gen_trihedron.ini' );
            [R_w_Cam_Trihedron, R_w_LRF_Trihedron, t_w_Rig_Trihedron, ~, ~] = generate_random_poses( Nobs, gen_config_file, Rig );
        end
        
        if WITHCORNER
            clearvars cornerOptim;
            cornerOptim = CCornerOptimization( K,...
                debug_level, maxIters,...
                minParamChange, minErrorChange);
            gen_config_file = fullfile( pwd, 'pose_gen_corner.ini' );
            [R_w_Cam_Corner, R_w_LRF_Corner, t_w_Rig_Corner, ~, ~] = generate_random_poses( Nobs, gen_config_file, Rig );
        end
        
        if WITHZHANG
            clearvars checkerOptim;
            checkerOptim = CCheckerboardOptimization( K,...
                debug_level, maxIters,...
                minParamChange, minErrorChange);
            gen_config_file = fullfile( pwd, 'pose_gen_checkerboard.ini' );
            [R_w_Cam_Checkerboard, R_w_LRF_Checkerboard, t_w_Rig_Checkerboard, ~, ~] = generate_random_poses( Nobs, gen_config_file, Rig );
        end

        % Store the observations
        for idx_N_co=1:Nobs
            % Correspondences for trihedron
                if WITHTRIHEDRON
                    % Update reference (Camera) pose in Rig for Trihedron
                    Rig.updateCamPose( R_w_Cam_Trihedron{idx_N_co}, t_w_Rig_Trihedron{idx_N_co} );
                    co = trihedron.getCorrespondence( Rig );
                    triOptim.stackObservation( co );
                    clear co
                end    
            
            % Correspondences for Kwak's algorithm
                if WITHCORNER
                    % Update reference (LRF) pose in Rig for Corner
                    Rig.updateLRFPose( R_w_LRF_Corner{idx_N_co}, t_w_Rig_Corner{idx_N_co} );
                    co = corner.getCorrespondence(Rig);
                    cornerOptim.stackObservation( co );
                    clear co
                end

                % Correspondences for Vasconcelos and Zhang's algorithm
                if WITHZHANG
                    % Update reference (LRF) pose in Rig for Checkerboard
                    Rig.updateLRFPose( R_w_LRF_Checkerboard{idx_N_co}, t_w_Rig_Checkerboard{idx_N_co} );
                    co = checkerboard.getCorrespondence( Rig );
                    checkerOptim.stackObservation( co );
                    clear co
                end
        end
        
        % Optimize only for value of Nsamples in N_co_n
        % --------------------------------------------------
        R0 = Rig.R_c_s;
        t0 = Rig.t_c_s;
        Rt0 = [R0 t0];
        for idx_N_co = 1:length(N_co_n)
            comp.setIndexes( idx_cam_sd, idx_scan_sd, idx_N_co, idx_Nsim );
            % Trihedron optimization
            if WITHTRIHEDRON
                triOptim.setNobs( N_co_n(idx_N_co) );
                triOptim.setInitialRotation( R0 );
                triOptim.setInitialTranslation( t0 );
                comp.TrihedronComp.optim( triOptim, WITHRANSAC );
            end
            
            % Kwak optimization
            if WITHCORNER                
                cornerOptim.setNobs( N_co_n(idx_N_co) );
                cornerOptim.setInitialRotation( R0 );
                cornerOptim.setInitialTranslation( t0 );
                comp.KwakComp.optim( cornerOptim );
            end
            
            % Zhang and Vasconcelos optimization
            if WITHZHANG
                checkerOptim.setNobs( N_co_n(idx_N_co) );
                checkerOptim.setInitialRotation( R0 );
                checkerOptim.setInitialTranslation( t0 );
                comp.VasconcelosComp.optim( checkerOptim );
                comp.ZhangComp.optim( checkerOptim );
            end            
        end
    end
end
end

toc
save(storeFile);