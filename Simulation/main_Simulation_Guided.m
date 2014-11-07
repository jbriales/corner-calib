% Main script for the synthetic simulations

% Simulated data for the Extrinsic Calibration of a 2D Lidar and a
% Monocular Camera based on Corner Structures without Pattern

% clear classes
clear;

% Main options:
main_sim_file = fullfile( pwd, 'main.ini' );
mainOpts = readConfigFile( main_sim_file, '[Simulation]' );
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
trihedron = CTrihedron( LTrihedron, eye(3), 0*[0 0 1]' );
corner = CCorner( LCorner, expmap( [-1 +1 0], deg2rad(-45) ) );
checkerboard = CCheckerboard( LCheckerboard, RotationZ(deg2rad(45))*RotationY(deg2rad(45)) );
pattern = { trihedron, corner, checkerboard };

% Generate Rig (Camera) poses for different patterns
gen_config_file = fullfile( pwd, 'pose_gen.ini' );
% Trihedron
[R_w_Cam_Trihedron, R_w_LRF_Trihedron, t_w_Rig_Trihedron, ~, ~] = ...
    generate_random_poses( Nobs, gen_config_file, '[Trihedron]', Rig );
% Corner
[R_w_Cam_Corner, R_w_LRF_Corner, t_w_Rig_Corner, ~, ~] =...
    generate_random_poses( Nobs, gen_config_file, '[Corner]', Rig );
% Checkerboard
[R_w_Cam_Checkerboard, R_w_LRF_Checkerboard, t_w_Rig_Checkerboard, ~, ~] =...
    generate_random_poses( Nobs, gen_config_file, '[Checkerboard]', Rig );

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
checkerOptim = CCheckerboardOptimization( K,...
    debug_level, maxIters,...
    minParamChange, minErrorChange);

for i=1:Nobs
    % Correspondences for Kwak's algorithm
    if WITHCORNER
        % Update reference (LRF) pose in Rig for Corner
        Rig.updateLRFPose( R_w_LRF_Corner{i}, t_w_Rig_Corner{i} );
        corr_ = corner.getCorrespondence(Rig);
        cornerOptim.stackObservation( corr_ );
    end
    
    % Correspondences for Vasconcelos and Zhang's algorithm
    % Update reference (LRF) pose in Rig for Checkerboard
    Rig.updateLRFPose( R_w_LRF_Checkerboard{i}, t_w_Rig_Checkerboard{i} );
    co = checkerboard.getCorrespondence( Rig );
    checkerOptim.stackObservation( co );
    
    % Correspondences for trihedron
    % Update reference (Camera) pose in Rig for Trihedron
    Rig.updateCamPose( R_w_Cam_Trihedron{i}, t_w_Rig_Trihedron{i} );
    co_ = trihedron.getCorrespondence( Rig );
    triOptim.stackObservation( co_ );
    
    if 0
        % Test multi-line framework
        Rwc = R_w_Cam_Trihedron{1};
        R0  = Rwc';
        to = {};
        Nbp = [];
        A_Nbp = {};
        for k=1:15
            Rig.updateCamPose( Rwc, t_w_Rig_Trihedron{k} );
            to_ = trihedron.getCorrespondence( Rig );
            if ~isempty(to_)
                to{end+1} = to_;
                Nbp = [Nbp to_.cam_reprN];
                A_Nbp{end+1} = to_.cam_A_reprN;
            end
        end
        ob = Manifold.S2( Nbp );
        obj_Nbp = Manifold.Dyn( ob );
        obj_Nbp.setRepresentationCov( blkdiag(A_Nbp{:}) );
        labels = kron(ones(1,obj_Nbp.Nvars/3),[1 2 3]); % For complete TO by now, but more flexible for future
        obj_Rtri = optimizeTrihedronNormals(obj_Nbp,labels,R0);
    end
    
    if WITHPLOTSCENE
        % Need to update Rig poses for plotting
        figure
        subplot(131)
        Rig.updateCamPose( R_w_Cam_Trihedron{i}, t_w_Rig_Trihedron{i} );
        trihedron.plotScene(Rig.Camera, Rig.Lidar);
        subplot(132)
        Rig.updateLRFPose( R_w_LRF_Corner{i}, t_w_Rig_Corner{i} );
        corner.plotScene(Rig.Camera, Rig.Lidar);
        subplot(133)
        Rig.updateLRFPose( R_w_LRF_Checkerboard{i}, t_w_Rig_Checkerboard{i} );
        checkerboard.plotScene(Rig.Camera, Rig.Lidar);
        set(gcf,'units','normalized','position',[0 0 1 1]);
        keyboard
        close
    end
end

% ------------- Trihedron ----------------
% Set number of observations to use
triOptim.setNobs(Nobs);
if WITHTRIHEDRON
    triOptim.setInitialRotation( [ 0 -1  0
                                   0  0 -1
                                   1  0  0 ] ); % Updated in RANSAC
    if WITHRANSAC
        triOptim.filterRotationRANSAC;
    end
    if WITHVERBOSE
        triOptim.disp_N_R_inliers;
    end
    R_c_s_w  = triOptim.optimizeRotation_Weighted;
    A_R_ort  = triOptim.FCov_R_W;

    if WITHRANSAC
        triOptim.filterTranslationRANSAC( Rig.R_c_s ); % Should receive some estimated rotation
    end
    if WITHVERBOSE
    triOptim.disp_N_t_inliers;
    end
    R0_for_t = R_c_s_w;
    triOptim.setInitialTranslation( Rig.t_c_s + 0.05*randn(3,1) );
    t_3D_w  = triOptim.optimizeTranslation_3D_Weighted( R0_for_t );
    A_t_3D = triOptim.FCov_t_3D_W;
    
    [R_global, t_global] = triOptim.optimizeGlobal_Ort_3D( R_c_s_w, t_3D_w );
    [R_global_3D, t_global_3D] = triOptim.optimizeGlobal_3D( R_c_s_w, t_3D_w );
    
    % Compute and check covariances
    if 0 % Monte Carlo for rotation
        keyboard
        % Create manifold input
        N = triOptim.cam_N;
        V = triOptim.LRF_V;
        cN = num2cell(N,1);
        cV = num2cell(V,1);
        for k=1:size(cN,2)
            obN{k} = Manifold.S2(cN{k});
            obV{k} = Manifold.S1(cV{k});
        end
        obData = Manifold.Dyn( obN{:}, obV{:} );
        % TODO yet
        out = Manifold.MonteCarloSim( ...,
            @(in)obj.simulateTOimage( reshape(in.X,2,4), Rig.Camera.sd ),...
            obj_pts, 'Ref', obj_xi, 'N', 1e4 );
        keyboard
    end
end

% ------------- Kwak -------------------
if WITHCORNER
    cornerOptim.setNobs(Nobs);
    if WITHVERBOSE
        cornerOptim.disp_N_obs;
    end
    % Generate random input (near GT)
    R_aux = Rig.R_c_s + randn(3,3)*0.08;
    [U,S,V] = svd(R_aux);
    Rt0 = [ U*V' , Rig.t_c_s + 0.05*randn(3,1) ];
    % Optimize
    cornerOptim.setInitialRotation( Rt0(1:3,1:3) );
    cornerOptim.setInitialTranslation( Rt0(1:3,4) );
    [R_k_cw, t_k_cw] = cornerOptim.optimizeRt_ConstWeighted;
end

% % ---------- Vasconcelos -------------------------
if WITHZHANG
    checkerOptim.setInitialRotation( Rig.R_c_s );
    checkerOptim.setInitialTranslation( Rig.t_c_s );
    [R_v,t_v] = checkerOptim.optimizeRt_Vasc;
    [R_z,t_z] = checkerOptim.optimizeRt_Zhang;
end

%% Compute hessian in convergence points for different methods
if WITHTRIHEDRON
    H_R = triOptim.FHes_Orthogonality( R_c_s_w );
    H_t_3D = triOptim.FHes_3D_PlaneDistance( Rig.R_c_s, t_3D_w );
    H_Rt = triOptim.FHes_Global_Ort_3D( R_global, t_global );
end
if WITHCORNER
    H_Rt_k = cornerOptim.FHes_2D_LineDistance( [R_k_cw, t_k_cw] );
end
if WITHPLOTHESSIAN
    figure
    f1 = 2; f2 = 4; index = 1;
    if WITHTRIHEDRON
        subplot(f1,f2,index), hold on, index = index + 1;
        eigv = eig(H_R); m = min(eigv);
        title(sprintf('Tri:R:W\n%i:%i:%i',round(eigv(1)/m),round(eigv(2)/m),round(eigv(3)/m)));
        plotcov3( zeros(3,1), H_R ); shading interp;
        subplot(f1,f2,index), hold on, index = index + 1;
        eigv = eig(H_t_3D); m = min(eigv);
        title(sprintf('Tri:t:3D:W\n%i:%i:%i',round(eigv(1)/m),round(eigv(2)/m),round(eigv(3)/m)));
        plotcov3( zeros(3,1), H_t_3D ); shading interp;
        subplot(f1,f2,index), hold on, index = index + 1;
        eigv = eig(H_t_2D); m = min(eigv);
        title(sprintf('Tri:t:2D:W\n%i:%i:%i',round(eigv(1)/m),round(eigv(2)/m),round(eigv(3)/m)));
        plotcov3( zeros(3,1), H_t_2D ); shading interp;
        subplot(f1,f2,index), hold on, index = index + 1;
        eigv = eig(H_Rt(1:3,1:3)); m = min(eigv);
        title(sprintf('Tri:Rt:Global:W\n%i:%i:%i',round(eigv(1)/m),round(eigv(2)/m),round(eigv(3)/m)));
        plotcov3( zeros(3,1), H_Rt(1:3,1:3) ); shading interp;
    end
    if WITHCORNER
        index = f2 + 1; % New row
        subplot(f1,f2,index), hold on, index = index + 1;
        eigv = eig(H_Rt_k(1:3,1:3)); m = min(eigv);
        title(sprintf('Kwak:R:NW\n%i:%i:%i',round(eigv(1)/m),round(eigv(2)/m),round(eigv(3)/m)));
        plotcov3( zeros(3,1), H_Rt_k(1:3,1:3) ); shading interp;
        subplot(f1,f2,index), hold on, index = index + 1;
        eigv = eig(H_Rt_k(4:6,4:6)); m = min(eigv);
        title(sprintf('Kwak:t:NW\n%i:%i:%i',round(eigv(1)/m),round(eigv(2)/m),round(eigv(3)/m)));
        plotcov3( zeros(3,1), H_Rt_k(4:6,4:6) ); shading interp;
        subplot(f1,f2,index), hold on, index = index + 1;
        eigv = eig(H_Rt_kC(1:3,1:3)); m = min(eigv);
        title(sprintf('Kwak:R:C:NW\n%i:%i:%i',round(eigv(1)/m),round(eigv(2)/m),round(eigv(3)/m)));
        plotcov3( zeros(3,1), H_Rt_kC(1:3,1:3) ); shading interp;
        subplot(f1,f2,index), hold on, index = index + 1;
        eigv = eig(H_Rt_kC(4:6,4:6)); m = min(eigv);
        title(sprintf('Kwak:t:C:NW\n%i:%i:%i',round(eigv(1)/m),round(eigv(2)/m),round(eigv(3)/m)));
        plotcov3( zeros(3,1), H_Rt_kC(4:6,4:6) ); shading interp;
    end
end

%% Plot cost functions near GT
if WITHPLOTCOST
    figure('Name','Trihedron Rotation: Orthogonality cost function');
    title('Trihedron Rotation: Orthogonality cost function');
    triOptim.plotRotationCostFunction( Rig.R_c_s );
    
    figure('Name','Trihedron Rotation: Global cost function');
    title('Trihedron Rotation: Global cost function');
    triOptim.plotRotation_Global_CostFunction( Rig.R_c_s, Rig.t_c_s );
    
    if WITHCORNER
        figure('Name','Corner Rotation: 2D distance cost function');
        title('Corner Rotation: 2D distance cost function');
        cornerOptim.plotRotationCostFunction( Rig.R_c_s, Rig.t_c_s );
        
        figure('Name','Corner Rotation: 2D distance cost function (only center)');
        title('Corner Rotation: 2D distance cost function (only center)');
        cornerOptim.plotRotation_C_CostFunction( Rig.R_c_s, Rig.t_c_s );
    end
    
    figure('Name','Trihedron Translation: 3D distance cost function');
    title('Trihedron Translation: 3D distance cost function');
    triOptim.plotTranslation_3D_CostFunction( Rig.R_c_s, Rig.t_c_s );
    
    figure('Name','Trihedron Translation: 2D distance cost function');
    title('Trihedron Translation: 2D distance cost function');
    triOptim.plotTranslation_2D_CostFunction( Rig.R_c_s, Rig.t_c_s );
    
    figure('Name','Trihedron Translation: Global cost function');
    title('Trihedron Translation: Global cost function');
    triOptim.plotTranslation_Global_CostFunction( Rig.R_c_s, Rig.t_c_s );
    
    if WITHCORNER
        figure('Name','Corner Translation: 2D distance cost function');
        title('Corner Translation: 2D distance cost function');
        cornerOptim.plotTranslationCostFunction( Rig.R_c_s, Rig.t_c_s );
        
        figure('Name','Corner Translation: 2D distance cost function (only center)');
        title('Corner Translation: 2D distance cost function (only center)');
        cornerOptim.plotTranslation_C_CostFunction( Rig.R_c_s, Rig.t_c_s );
    end
end

%% ---------- Display the errors -------------------------
if WITHVERBOSE
    fprintf('(*) -> best current method\n');
    
    if WITHTRIHEDRON
        fprintf('=============================================================\n');
        fprintf('(*) Trihedron (weighted) rotation error (deg): \t \t %f \n',...
            angularDistance(R_c_s_w,Rig.R_c_s) );
        fprintf('Trihedron (global W, 3D) rotation error (deg): \t \t %f \n',...
            angularDistance(R_global,Rig.R_c_s) );
        fprintf('Trihedron (global W, 3D only) rotation error (deg): \t %f \n',...
            angularDistance(R_global_3D,Rig.R_c_s) );
        fprintf('=============================================================\n');
        fprintf('Trihedron (    weighted, 3D) translation error (cm): \t %f \n',...
            norm(t_3D_w-Rig.t_c_s)*100 );
        fprintf('Trihedron (    global W, 3D) translation error (cm): \t %f \n',...
            norm(t_global-Rig.t_c_s)*100 );
        fprintf('Trihedron (    global W, 3D only) translation error (cm):%f \n',...
            norm(t_global_3D-Rig.t_c_s)*100 );
    end
    
    if WITHCORNER
        fprintf('=============================================================\n');
        fprintf('=============================================================\n');
        fprintf('(*) Kwak (const-weighted) rotation error (deg): \t %f \n',   angularDistance(R_k_cw,Rig.R_c_s ));
        fprintf('=============================================================\n');
        fprintf('(*) Kwak (const-weighted) translation error (cm): \t %f \n', 100*norm(t_k_cw - Rig.t_c_s) );
    end
    
    if WITHZHANG
        fprintf('=============================================================\n');
        fprintf('=============================================================\n');
        fprintf('Vasconcelos rotation error (deg): \t\t\t %f \n', angularDistance(R_v,Rig.R_c_s) );
        fprintf('Zhang rotation error (deg): \t\t\t\t %f \n', angularDistance(R_z,Rig.R_c_s) );
        fprintf('=============================================================\n');
        fprintf('Vasconcelos translation error (cm): \t\t\t %f \n', 100 * norm(t_v - Rig.t_c_s) );
        fprintf('Zhang translation error (cm): \t\t\t\t %f \n', 100 * norm(t_z - Rig.t_c_s) );
    end
    toc
end