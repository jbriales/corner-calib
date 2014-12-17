close all

% WITH_ELQURSH = false;
WITH_ELQURSH = true;
WITH_TRIHEDRON = true;

% WITH_PLOT = false;
WITH_PLOT = true;

% Read camera  options
rig_config_file = fullfile( pwd, 'rig.ini' );
rigOpts = readConfigFile( rig_config_file );
extractStructFields( rigOpts );

% Generate random pose for camera
gen_conf_F = fullfile( pwd, 'pose_gen.ini' );
poseFactory = CRandomPoses( readConfigFile( gen_conf_F,'[Vanishing]' ) );
[R_w_c, t_w_c] = poseFactory.gen( 1 );

% Create simulated camera
% TEMPORAL: cam_sd
cam_sd = 0;
Cam = CSimCamera( R_w_c, t_w_c, K, res, f, cam_sd );

% Create pattern
Cube_sd = 0.001;
Cube = CHalfCube( 1, eye(3), zeros(3,1), Cube_sd );

% Plot scene
if WITH_PLOT
hF_scene = figure;
Cube.plotScene(Cam);
end

% Compute projection of points and lines
% [uv_proj, uv_pixels] = Cube.getProjection( Cam );
% lines = Cube.getProjectionLines( Cam );

R_gt = R_w_c';

if WITH_ELQURSH
% Solve with Elqursh
primitives = Cube.elqursh_primitives;
Nprim = length(primitives);
K = Cam.K;
err_elqursh = zeros(1,length(primitives));
I = eye(3);
cell_d = cell(1,Nprim);
cell_v = cell(1,Nprim);
hF_image = figure; hold on;
for k=1:Nprim
    prim = primitives{k};
    
    if WITH_PLOT
    figure(hF_scene);
    h_prim = Cube.plot_prim( prim );
    end

    prim2D = prim.project(Cam);
    if WITH_PLOT
    figure(hF_image);
    prim2D.plot;
    end
    
    for ii=1:3
        line{ii} = prim(ii).projectLine( Cam );
    end
    v2 = cross( line{2}, line{3} );
    v1 = null( [ v2'*(inv(K)'*inv(K)) ; line{1}' ] );
    
%     dir1 = find(prim(1).v~=0);
%     dir2 = find(prim(2).v~=0);
%     dir3 = setdiff(1:3,[dir1 dir2]);
    [~,dir1] = max(abs(prim(1).v));
    [~,dir2] = max(abs(prim(2).v));
    dir3 = setdiff(1:3,[dir1 dir2]);    

    R = eye(3);
    R(:,dir1) = snormalize(K \ v1);
    R(:,dir2) = snormalize(K \ v2);
    switch dir3
        case 1
            R(:,1) = cross(R(:,2),R(:,3));
        case 2
            R(:,2) = cross(R(:,3),R(:,1));
        case 3
            R(:,3) = cross(R(:,1),R(:,2));
    end
    
    % Compute 4 possible sign combinations
    CR = cell(1,4);
    CR{1} = R;
    CR{2} = R * diag([-1 1 -1]);
    CR{3} = R * diag([-1 -1 1]);
    CR{4} = R * diag([1 -1 -1]);
    % Find good one minimizing Frobenius norm
    d = zeros(1,4);
    for iR = 1:4
        d(iR) = angularDistance(CR{iR},R_gt);
    end
    [~,iR] = min(d);
    R = CR{iR};
    % Error
    err_elqursh(k) = angularDistance(R,R_gt);
    
    % Store vanishing points and directions for global Procrustes
    cell_d{k} = I(:,[dir1 dir2]);
%     cell_d{k} = [prim(1).v, prim(2).v];
    cell_v{k} = R(:,[dir1 dir2]); % Equivalent to inv(K)*K*rk
end
% Set camera image borders
if WITH_PLOT
figure(hF_image);
Cam.setImageBorder;
end

D1 = [cell_d{:}];
D2 = [cell_v{:}];
norm(D1-R'*D2,'fro')
% Compute through Procrustes problem:
[U,S,V] = svd( D2*D1' );
R_pro = U*V';
fprintf('Elqursh Procrustes error: %f\n', angularDistance(R_pro,R_gt));

figure
hist( err_elqursh );

% keyboard
end


if WITH_TRIHEDRON
% Compute with Trihedron method
primitives = Cube.trihedron_primitives;
Nprim = length(primitives);

K = Cam.K;
err_trihedron = cell(1,length(primitives));
cell_d = cell(1,Nprim);
cell_v = cell(1,Nprim);
if WITH_PLOT
hF_image = figure; hold on;
end
for k=1:Nprim
% for k=[4 14 18 20]
% for k=14
% for k=19:21
    prim = primitives{k};
    % Take care of segments orientation
    % (building a dextrorotatory system)
    for kk=1:3
        % Valid only for canonical directions (simulation)
        if any(prim(kk).v < 0)
            prim(kk) = prim(kk).inverse;
        end
    end
    
    if WITH_PLOT
    figure(hF_scene);
    h_prim = Cube.plot_prim( prim );
    end

    prim2D = prim.project(Cam);
    if WITH_PLOT
    figure(hF_image);
    prim2D.plot('-k',{'x','y','z'});
    end
    
    for ii=1:3
        C_Nbp{ii} = snormalize( K' * prim(ii).projectLine( Cam ) );
    end
    Nbp = [C_Nbp{:}];
    
    % Current solution only works for planes intersecting in one line!
    % Generic position not valid (lines intersecting in one single point)
%     disp( det( Nbp ) );
%     if rank(Nbp) == 2 || abs(det(Nbp)) < 1e-2
        trihedronSolver = CTrihedronSolver( Nbp, K );
        trihedronSolver.loadSegments( prim2D );
        V_tri = trihedronSolver.solve;
        if WITH_PLOT
            plotHomLineWin(inv(K')*trihedronSolver.Nbp(:,3), 'r');
        end
%         V_tri = trihedronSolver.solve_withGT( R_gt );
        
        % Error
        err_trihedron{k} = angularDistance(V_tri,R_gt);
        
        % Store vanishing points and directions for global Procrustes
        cell_d{k} = eye(3);
        cell_v{k} = V_tri; % Equivalent to inv(K)*K*rk
%     end
end
% Set camera image borders
if WITH_PLOT
figure(hF_image);
Cam.setImageBorder;
end

D1 = [cell_d{:}];
D2 = [cell_v{:}];
% norm(D1-R'*D2,'fro')
% Compute through Procrustes problem:
[U,S,V] = svd( D2*D1' );
R_pro = U*V';
fprintf('Trihedron Procrustes error: %f\n', angularDistance(R_pro,R_gt));

figure
hist( [err_trihedron{:}] );

end

if WITH_TRIHEDRON && WITH_ELQURSH
    figure; title('Error histogram');hold on;
    plot(1, err_elqursh, '*b');
    plot(2, [err_trihedron{:}], '*g');
    ax = axis; axis([0 3 ax(3:4)]);
    
    figure; title('Error histogram');hold on;
    hist( err_elqursh, 20, 0:0.1:10 );
    hist( [err_trihedron{:}], 0:0.1:10 );
    h = findobj(gca,'Type','patch');
    set(h(2),'FaceColor','b','EdgeColor','k');
    set(h(1),'FaceColor','g','EdgeColor','k');
end