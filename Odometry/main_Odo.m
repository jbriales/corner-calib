% Read camera options
rig_config_file = fullfile( pwd, 'rig.ini' );
rigOpts = readConfigFile( rig_config_file );
extractStructFields( rigOpts );
clear rigOpts

% Generate random pose for camera
gen_conf_F = fullfile( pwd, 'pose_gen.ini' );
[R_w_c, ~, t_w_c, ~, ~] = ...
                generate_random_poses( 1, gen_conf_F, '[Trihedron]' );

% Create simulated camera
Cam = CSimCamera( R_w_c, t_w_c, K, res, f, cam_sd );

% Create pattern
cube = CCube( 1, eye(3), zeros(3,1) );

% Plot scene
% cube.plotScene(Cam);

% Compute projection of points and lines
% [uv_proj, uv_pixels] = cube.getProjection( Cam );
lines = cube.getProjectionLines( Cam );
keyboard % For final check

% Solve with Elqursh
v2 = cross( lines.x{1}, lines.x{2} );
v1 = null( [ v2' * (inv(K')*inv(K)); lines.z{1}' ] );
r1 = K \ v1;
r2 = K \ v2;
r3 = cross(r1,r2);
R  = [r1 r2 r3];