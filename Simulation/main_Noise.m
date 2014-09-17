function  comp = main_Noise( )
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
% noise_config_file = fullfile( pwd, 'sim_noise.ini' );
noiseOpts = readConfigFile( noise_config_file );
extractStructFields( noiseOpts );
clear noiseOpts

cam_sd_n  =   logspace(cam_sd_range(1),cam_sd_range(2),cam_sd_N);
scan_sd_n = 3*logspace(scan_sd_range(1),scan_sd_range(2),scan_sd_N);
N_co_N    = size(N_co_n,2);

% Create the output structure for comparison
comp = CComparison( Rig.Rt_c_s, Nsim, cam_sd_n, scan_sd_n, N_co_n );

% Create objects for storage and optimization
optim_config_file = fullfile( pwd, 'optim_config.ini' );
optimOpts = readConfigFile( optim_config_file );
extractStructFields( optimOpts );
clear optimOpts

% Start the loop
tic
Nobs = N_co_n(end); %#ok<COLND>
counter_loops = 1;
total_loops   = numel(cam_sd_n) * numel(scan_sd_n) * numel(N_co_n) * Nsim;
first_time    = toc;
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
            [R_w_Cam_Trihedron, ~, t_w_Rig_Trihedron, ~, ~] = generate_random_poses( Nobs, gen_config_file, Rig );
        end
        
        if WITHCORNER
            clearvars cornerOptim;
            cornerOptim = CCornerOptimization( K,...
                debug_level, maxIters,...
                minParamChange, minErrorChange);
            gen_config_file = fullfile( pwd, 'pose_gen_corner.ini' );
            [~, R_w_LRF_Corner, t_w_Rig_Corner, ~, ~] = generate_random_poses( Nobs, gen_config_file, Rig );
        end
        
        if WITHZHANG
            clearvars checkerOptim;
            checkerOptim = CCheckerboardOptimization( K,...
                debug_level, maxIters,...
                minParamChange, minErrorChange);
            gen_config_file = fullfile( pwd, 'pose_gen_checkerboard.ini' );
            [~, R_w_LRF_Checkerboard, t_w_Rig_Checkerboard, ~, ~] = generate_random_poses( Nobs, gen_config_file, Rig );
        end

        % Store the observations
        for idx_N_co=1:Nobs
            % Correspondences for trihedron
                if WITHTRIHEDRON
                    % Update reference (Camera) pose in Rig for Trihedron
                    Rig.updateCamPose( R_w_Cam_Trihedron{idx_N_co}, t_w_Rig_Trihedron{idx_N_co} );
                    co = trihedron.getCorrespondence( Rig );
                    if isempty(co)
                        co = getCorrectCorrespondence( trihedron, Rig );
                    end
                    triOptim.stackObservation( co );
                    clear co
                end    
            
            % Correspondences for Kwak's algorithm
                if WITHCORNER
                    % Update reference (LRF) pose in Rig for Corner
                    Rig.updateLRFPose( R_w_LRF_Corner{idx_N_co}, t_w_Rig_Corner{idx_N_co} );
                    co = corner.getCorrespondence(Rig);
                    if isempty(co)
                        co = getCorrectCorrespondence( corner, Rig );
                    end
                    cornerOptim.stackObservation( co );
                    clear co
                end

                % Correspondences for Vasconcelos and Zhang's algorithm
                if WITHZHANG
                    % Update reference (LRF) pose in Rig for Checkerboard
                    Rig.updateLRFPose( R_w_LRF_Checkerboard{idx_N_co}, t_w_Rig_Checkerboard{idx_N_co} );
                    co = checkerboard.getCorrespondence( Rig );
                    if isempty(co)
                        co = getCorrectCorrespondence( checkerboard, Rig );
                    end
                    checkerOptim.stackObservation( co );
                    clear co
                end
                
                if WITHPLOTSCENE
                    % Need to update Rig poses for plotting
                    figure
                    if WITHTRIHEDRON
                        subplot(131)
                        Rig.updateCamPose( R_w_Cam_Trihedron{idx_N_co}, t_w_Rig_Trihedron{idx_N_co} );
                        trihedron.plotScene(Rig.Camera, Rig.Lidar);
                    end
                    if WITHCORNER
                        subplot(132)
                        Rig.updateLRFPose( R_w_LRF_Corner{idx_N_co}, t_w_Rig_Corner{idx_N_co} );
                        corner.plotScene(Rig.Camera, Rig.Lidar);
                    end
                    if WITHZHANG
                        subplot(133)
                        Rig.updateLRFPose( R_w_LRF_Checkerboard{idx_N_co}, t_w_Rig_Checkerboard{idx_N_co} );
                        checkerboard.plotScene(Rig.Camera, Rig.Lidar);
                    end
                    set(gcf,'units','normalized','position',[0 0 1 1]);
                    keyboard
                    close
                end
        end
        
        % Optimize only for value of Nsamples in N_co_n
        % --------------------------------------------------
        R0 = Rig.R_c_s;
        t0 = Rig.t_c_s;
        Rt0 = [R0 t0];
        for idx_N_co = 1:length(N_co_n)
            fprintf('Optimization %d of %d\n',counter_loops,total_loops);
            counter_loops = counter_loops + 1;
            
            % Set static indexes common for all comp objects
            CBaseComp.idx_cam_sd( idx_cam_sd );
            CBaseComp.idx_scan_sd( idx_scan_sd );
            CBaseComp.idx_N_co( idx_N_co );
            CBaseComp.idx_Nsim( idx_Nsim );
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
            
            % Estimate duration of simulation
%             mean_time = (toc-first_time)/counter_loops;
            mean_time = (toc-first_time)/counter_loops;
            fprintf('Mean time: %.2f\tEstimated left duration: %.2f\n',mean_time,(total_loops-counter_loops)*mean_time);
            first_time = toc;
        end
    end
end
    
end
save(storeFile); % For further repetition

% Plot results in boxplots
comp.plotCameraNoise;
end

function co = getCorrectCorrespondence( pattern, Rig )
co = [];
counter = 0;

if isa(pattern,'CTrihedron')
    gen_config_file = fullfile( pwd, 'pose_gen_trihedron.ini' );
elseif isa(pattern,'CCorner')
    gen_config_file = fullfile( pwd, 'pose_gen_corner.ini' );
elseif isa(pattern,'CCheckerboard')
    gen_config_file = fullfile( pwd, 'pose_gen_checkerboard.ini' );
end

while isempty(co) && counter < 500

[R_w_Cam, R_w_LRF, t_w_Rig, ~, ~] = generate_random_poses( 1, gen_config_file, Rig );

% Update reference (LRF or Cam) pose in Rig for Pattern
if isa(pattern,'CTrihedron')
    Rig.updateCamPose( R_w_Cam, t_w_Rig );
elseif isa(pattern,'CCorner')
    Rig.updateLRFPose( R_w_LRF, t_w_Rig );
elseif isa(pattern,'CCheckerboard')
    Rig.updateLRFPose( R_w_LRF, t_w_Rig );
end

co = pattern.getCorrespondence(Rig);
end
if counter >= 500
    warning('Could not generate a correct observation');
else
    cprintf('green', 'Solved empty observation\n')
end
end