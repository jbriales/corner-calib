% Main script for the synthetic simulations

% Simulated data for the Extrinsic Calibration of a 2D Lidar and a
% Monocular Camera based on Corner Structures without Pattern

clear classes

% Generate Rig (Camera) poses
[R_w_c, t_w_c] = generate_random_poses( );
% N = length(R_w_c);
Nsamples = 1000;

% Set Rig properties
R_c_s = [ 0 -1 0 ; 0 0 -1 ; 1 0 0 ] * RotationZ(deg2rad(30));
t_c_s = [0.5 0.25 0]';
N = 1081; FOVd = 270.2; scan_sd = 0.03; d_range = [0.1 30];
K = [ 1050 0 480
      0 1050 270
      0    0   1 ];
res = [960 540]; f = 1; cam_sd = 1;
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
        trihedron.plotScene( Rig.Lidar, Rig.Camera );
        corner.plotScene( Rig.Lidar, Rig.Camera );
        checkerboard.plotScene( Rig.Lidar, Rig.Camera );
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