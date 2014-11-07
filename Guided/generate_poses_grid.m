function [R, t] = generate_poses_grid( dx, dang, height, d_min, d_max, kFOV, Rig )
% [R, t] = generate_poses_grid( dx, dang, height, d_min, d_max, kFOV, Rig )

[X,Y,Z] = meshgrid( d_min:dx:d_max, d_min:dx:d_max, d_min:dx:height );

% World to Cam translation
t0 = [ X(:) Y(:) Z(:) ]';

% World to Cam rotation
% Set firstly camera pointing to World origin with X axis in XY world plane
% direction
R_z = -snormalize(t0);
R_x = snormalize( skew(-canvec(3)) * R_z );
R_y = cross(R_z,R_x,1);
R0 = reshape( [R_x;R_y;R_z], 3,3,[] );

% Set of possible rotation variations:
drad = deg2rad( dang );
span_z = -pi/2:drad:+pi/2;
span_y = -kFOV * Rig.Camera.FOVh : drad : +kFOV * Rig.Camera.FOVh;
span_x = -kFOV * Rig.Camera.FOVv : drad : +kFOV * Rig.Camera.FOVv;
[angX,angY,angZ] = meshgrid( span_x, span_y, span_z );
urot = [angX(:), angY(:), angZ(:)]';
inc_rot = expmap( urot );

% Get all possible subrotations for each position
R = zeros(3,3,0);
for k=1:size(R0,3)
    R_ = multiprod( R0(:,:,k), inc_rot, [1 2], [1 2] );
    R = cat( 3, R, R_ );
end
% Repeat translations for each possible rotation
t = repmat( t0, size(inc_rot,3), 1 );
t = reshape( t, 3,[] );

N = size(R,3);
% Transform output to cell arrays:
t = mat2cell( t, 3, ones(1,N) );
R = mat2cell( R, 3, 3, ones(1,N) );
R = permute( R, [1 3 2] );

end