% Main script for the synthetic simulations

% Simulated data for the Extrinsic Calibration of a 2D Lidar and a
% Monocular Camera based on Corner Structures without Pattern

% clear classes
clear;

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
    debug_level, maxIters,...
    minParamChange, minErrorChange);
cornerOptim = CCornerOptimization( K,...
    debug_level, maxIters,...
    minParamChange, minErrorChange);

for i=1:Nsamples
    % Update reference (Camera) pose in Rig
    Rig.updatePose( R_w_c{i}, t_w_c{i} );
        
    % Correspondences for Kwak's algorithm
    corr_ = corner.getCorrespondence(Rig); 
    cornerOptim.stackObservation( corr_ );
    for j = 1:2
        for k = 1:3
            corner_corresp{j,k,i} = corr_{j,k};
        end
    end
%     figure
%     corner.plotScene(Rig.Camera, Rig.Lidar);
%     set(gcf,'units','normalized','position',[0 0 1 1]);
%     close;
    
    % Correspondences for Vasconcelos and Zhang's algorithm
    check_corresp{1,i} = checkerboard.p2D; 
    check_corresp{2,i} = checkerboard.getProjection( Rig.Camera );    
    check_corresp{3,i} = 1000 * cell2mat(checkerboard.getScan( Rig.Lidar ));
    
    % Correspondences for trihedron
%     figure
%     trihedron.plotScene(Rig.Camera, Rig.Lidar);
%     set(gcf,'units','normalized','position',[0 0 1 1]);
%     close;
    co_ = trihedron.getCorrespondence( Rig );
%     co(i) = co_;
    triOptim.stackObservation( co_ );
end

% ------------- Trihedron ----------------
triOptim.setInitialRotation( [ 0 -1  0
                               0  0 -1
                               1  0  0 ] ); % Updated in RANSAC
triOptim.filterRotationRANSAC;
triOptim.disp_N_R_inliers;
R_c_s_nw = triOptim.optimizeRotation_NonWeighted;
R_c_s_dw = triOptim.optimizeRotation_DiagWeighted;
R_c_s_w  = triOptim.optimizeRotation_Weighted;

% Plot rotation cost function near GT
% triOptim.plotRotationCostFunction( Rig.R_c_s );
                           
% co0 = co;
% s_FinalOptimization
% R0 = R_c_s;
% solveTranslation
% solveTranslation_3D
R_c_s
R_c_s_w
R_c_s_dw
R_c_s_nw
angularDistance( R_c_s, R_c_s_w )
angularDistance( R_c_s, R_c_s_dw )

triOptim.filterTranslationRANSAC( Rig.R_c_s ); % Should receive some estimated rotation
triOptim.disp_N_t_inliers;
triOptim.setInitialTranslation( Rig.t_c_s + 0.1*randn(3,1) );
t_3D_nw = triOptim.optimizeTranslation_3D_NonWeighted( Rig.R_c_s );
t_3D_w  = triOptim.optimizeTranslation_3D_Weighted( Rig.R_c_s );
t_2D_nw = triOptim.optimizeTranslation_2D_NonWeighted( Rig.R_c_s );
t_3D_nw
t_3D_w
t_2D_nw

% Plot translation cost function near GT
% triOptim.plotTranslation_3D_CostFunction( Rig.R_c_s, Rig.t_c_s );


% ------------- Kwak -------------------
% Generate random input (near GT)
R_aux = Rig.R_c_s + randn(3,3)*0.01;
[U,S,V] = svd(R_aux);
Rt0 = [ U*V' , Rig.t_c_s + 0.05*randn(3,1) ];
% Optimize
Rt_knw  = corner.optim(corner_corresp,Rt0,0,Rig);
Rt_kw   = corner.optim(corner_corresp,Rt0,1,Rig);
Rt_gt = [Rig.R_c_s Rig.t_c_s];
% Check distance between optimization and initial point
fprintf('Kwak (NW): Change in rotation = %f\n',angularDistance(Rt_knw(1:3,1:3),Rt0(1:3,1:3)));
fprintf('Kwak ( W): Change in rotation = %f\n',angularDistance(Rt_kw(1:3,1:3),Rt0(1:3,1:3)));
fprintf('Kwak (NW): Change in translation = %f\n',norm(Rt_knw(1:3,1:3) - Rt0(1:3,1:3)));
fprintf('Kwak ( W): Change in translation = %f\n',norm(Rt_kw(1:3,1:3)  - Rt0(1:3,1:3)));
fprintf('\n\n')

% % ---------- Vasconcelos -------------------------
% [T_planes,lidar_points] = checkerboard.getCalibPlanes( Rig, check_corresp );
% [T, ~,~,~,~] = lccMinSol(T_planes,lidar_points);
% [T_z, ~,~,~,~] = lccZhang(T_planes, lidar_points);
% x_v = pose_inverse(T); x_v(1:3,4) = x_v(1:3,4)/1000;
% x_z = pose_inverse(T_z); x_z(1:3,4) = x_z(1:3,4)/1000;

% ---------- Display the errors -------------------------


%fprintf('Kwak translation error (m): \t %f \n', norm(x_w(:,4) - x_gt(:,4)) );
fprintf('Trihedron (weighted) rotation error (deg): \t \t %f \n',...
    angularDistance(R_c_s_w,Rig.R_c_s) );
fprintf('Trihedron (diag-weighted) rotation error (deg): \t %f \n',...
    angularDistance(R_c_s_dw,Rig.R_c_s) );
fprintf('Trihedron (non-weighted) rotation error (deg): \t \t %f \n',...
    angularDistance(R_c_s_nw,Rig.R_c_s) );
fprintf('Trihedron (non-weighted, 3D) translation error (cm): \t %f \n',...
    norm(t_3D_nw-Rig.t_c_s)*100 );
fprintf('Trihedron (    weighted, 3D) translation error (cm): \t %f \n',...
    norm(t_3D_w-Rig.t_c_s)*100 );
fprintf('Trihedron (non-weighted, 2D) translation error (cm): \t %f \n',...
    norm(t_3D_nw-Rig.t_c_s)*100 );

fprintf('Kwak (weighted) rotation error (deg): \t \t \t %f \n', angularDistance(Rt_kw(1:3,1:3),Rt_gt(1:3,1:3)) );
fprintf('Kwak (non-weighted) rotation error (deg): \t \t %f \n', angularDistance(Rt_knw(1:3,1:3),Rt_gt(1:3,1:3)) );
fprintf('Kwak (weighted) translation error (cm): \t\t %f \n', 100*norm(Rt_kw(:,4) - Rt_gt(:,4)) );
fprintf('Kwak (non-weighted) translation error (cm): \t\t %f \n', 100*norm(Rt_knw(:,4) - Rt_gt(:,4)) );

% fprintf('Vasconcelos translation error (cm): \t %f \n', 100 * norm(x_v(1:3,4) - x_gt(1:3,4)) );
% fprintf('Vasconcelos rotation error (deg): \t %f \n', angularDistance(x_v(1:3,1:3),x_gt(1:3,1:3)) );

% fprintf('Zhang translation error (cm): \t %f \n', 100 * norm(x_z(1:3,4) - x_gt(1:3,4)) );
% fprintf('Zhang rotation error (deg): \t %f \n', angularDistance(x_z(1:3,1:3),x_gt(1:3,1:3)) );

toc