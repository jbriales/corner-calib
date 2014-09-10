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
Nsamples = length(R_w_c);
corresp  = cell(2,3,Nsamples);

% Set Rig properties
% R_c_s = [ 0 -1 0 ; 0 0 -1 ; 1 0 0 ] * RotationZ(deg2rad(30));
% t_c_s = [0.5 0.25 0]';
R_c_s = [ 0 -1 0 ; 0 0 -1 ; 1 0 0 ];
t_c_s = [0.15 0 0]';

N = 1081; FOVd = 270.2; scan_sd = 0.000000000000000000003; d_range = [0.1 30];
K = [ 1050 0 480
      0 1050 270
      0    0   1 ];
res = [960 540]; f = 1; cam_sd = 0.0000000000000000000000001;
Rig = CSimRig( eye(3), zeros(3,1), R_c_s, t_c_s,... % Extrinsic options
               N, FOVd, scan_sd, d_range,... % Lidar options
               K, res, f, cam_sd ); % Camera options 
           
% Rig = CSimRig( eye(3), zeros(3,1), R_c_s, t_c_s,... % Extrinsic options
%    N, FOVd, 0, d_range,... % Lidar options
%    K, res, f, 0 ); % Camera options       

trihedron = CTrihedron;
corner = CCorner( expmap( [-1 +1 0], deg2rad(-45) ) );
checkerboard = CCheckerboard( RotationZ(deg2rad(45))*RotationY(deg2rad(45)) );
pattern = { trihedron, corner, checkerboard };

corner_corresp = cell(2,3,Nsamples);

tic
for i=1:Nsamples
    % Update reference (Camera) pose in Rig
    Rig.updatePose( R_w_c{i}, t_w_c{i} );
    
%     if 1 % For plotting
% %         trihedron.getCornerData( Rig.Camera );
% %         trihedron.plotScene( Rig.Lidar, Rig.Camera );
%         corner.plotScene( Rig.Lidar, Rig.Camera );
% %         checkerboard.plotScene( Rig.Lidar, Rig.Camera );
%         set(gcf,'units','normalized','position',[0 0 1 1]);
% %         pause( )
% %         close
%     else % For computation only
%         tic
%         for j=1:3
%             pattern{j}.getProjection( Rig.Camera );
%             pattern{j}.getScan( Rig.Lidar );
%         end
%         toc
%     end
    
    % Correspondences for Kwak's algorithm
    corr_ = corner.getCorrespondence(Rig); 
    for j = 1:2
        for k = 1:3
            corner_corresp{j,k,i} = corr_{j,k};
        end
    end   
    
    % Correspondences for Vasconcelos and Zhang's algorithm
    check_corresp{1,i} = checkerboard.p2D; 
    check_corresp{2,i} = checkerboard.getProjection( Rig.Camera );    
    check_corresp{3,i} = 1000 * cell2mat(checkerboard.getScan( Rig.Lidar ));
    
    % Correspondences for trihedron 
    co(i) = trihedron.getCorrespondence( Rig );
    
    
end

% ------------- Trihedron ----------------
co0 = co;
s_FinalOptimization
R0 = R_c_s;
% solveTranslation
% solveTranslation_3D
R_c_s
R_c_s_w
R_c_s_nw
angularDistance( R_c_s, R_c_s_w )

   
% ------------- Kwak -------------------
x_w  = corner.optim(corner_corresp,0,1,Rig);
x_gt = [Rig.R_c_s Rig.t_c_s];

% ---------- Vasconcelos -------------------------
[T_planes,lidar_points] = checkerboard.getCalibPlanes( Rig, check_corresp );
[T, ~,~,~,~] = lccMinSol(T_planes,lidar_points);
[T_z, ~,~,~,~] = lccZhang(T_planes, lidar_points);
x_v = pose_inverse(T); x_v(1:3,4) = x_v(1:3,4)/1000;
x_z = pose_inverse(T_z); x_z(1:3,4) = x_z(1:3,4)/1000;

% ---------- Display the errors -------------------------
%fprintf('Kwak translation error (m): \t %f \n', norm(x_w(:,4) - x_gt(:,4)) );
fprintf('Trihedron rotation error (deg): \t %f \n', angularDistance(R_c_s_w(1:3,1:3),x_gt(1:3,1:3)) );

fprintf('Kwak translation error (m): \t %f \n', norm(x_w(:,4) - x_gt(:,4)) );
fprintf('Kwak rotation error (deg): \t %f \n', angularDistance(x_w(1:3,1:3),x_gt(1:3,1:3)) );

fprintf('Vasconcelos translation error (m): \t %f \n', norm(x_v(1:3,4) - x_gt(1:3,4)) );
fprintf('Vasconcelos rotation error (deg): \t %f \n', angularDistance(x_v(1:3,1:3),x_gt(1:3,1:3)) );

fprintf('Zhang translation error (m): \t %f \n', norm(x_z(1:3,4) - x_gt(1:3,4)) );
fprintf('Zhang rotation error (deg): \t %f \n', angularDistance(x_z(1:3,1:3),x_gt(1:3,1:3)) );

toc