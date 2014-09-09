% Main script for the synthetic simulations

% Simulated data for the Extrinsic Calibration of a 2D Lidar and a
% Monocular Camera based on Corner Structures without Pattern

% clear classes

% Generate Rig (Camera) poses
% [R_w_c, t_w_c] = generate_random_poses( );
gen_config_file = fullfile( pwd, 'pose_gen.ini' );
[R_w_c, t_w_c, rand_ang_z, rand_ang_x] = generate_random_poses( gen_config_file ); % For debug purposes only
rand_ang_x = rad2deg( rand_ang_x );
rand_ang_z = rad2deg( rand_ang_z );
% N = length(R_w_c);
Nsamples = 1000;

% Set Rig properties
rig_config_file = fullfile( pwd, 'rig.ini' );
rigOpts = readConfigFile( rig_config_file );
extractStructFields( rigOpts );
clear rigOpts
Rig = CSimRig( eye(3), zeros(3,1), R_c_s, t_c_s,... % Extrinsic options
               N, FOVd, scan_sd, d_range,... % Lidar options
               K, res, f, cam_sd ); % Camera options                

trihedron = CTrihedron;
corner = CCorner( expmap( [-1 +1 0], deg2rad(-45) ) );
checkerboard = CCheckerboard( RotationZ(deg2rad(45))*RotationY(deg2rad(45)) );
pattern = { trihedron, corner, checkerboard };
tic
for i=1:Nsamples
    % Update reference (Camera) pose in Rig
    Rig.updatePose( R_w_c{i}, t_w_c{i} );
    
    if 1 % For plotting
%         trihedron.getCornerData( Rig.Camera );
%         trihedron.plotScene( Rig.Lidar, Rig.Camera );
        corner.plotScene( Rig.Lidar, Rig.Camera );
%         checkerboard.plotScene( Rig.Lidar, Rig.Camera );
        set(gcf,'units','normalized','position',[0 0 1 1]);
        pause( )
        close
    else % For computation only
        tic
        for j=1:3
            pattern{j}.getProjection( Rig.Camera );
            pattern{j}.getScan( Rig.Lidar );
        end
        toc
    end
end
toc