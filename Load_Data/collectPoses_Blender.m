function [ gt, v ] = collectPoses_Blender( gtPath )
%[ gt ] = collectPoses( gtPath )
%   Input:
%       gtPath - groundtruth.txt path
% 
%   Output: 
%       gt - STRUCT with groundtruth data
%          * ts - frame of the groundtruth   [1xN]
%          * t - displacement [tx,ty,tz]'       [3xN]
%          * R  - rotation matrix               [3x3xN]
%          * u  - angle-axis rotation           [3xN]


FID = fopen(fullfile(gtPath,'groundtruth.txt'));
if FID < 0
    error('EXISTENCE: Check that groundtruth.txt exists in dataset folder')
end
fgetl( FID ); fgetl( FID ); fgetl( FID ); % Read the first 3 rows

v = fscanf(FID, '%f %f %f %f %f %f %f\n', [7 Inf]);

fclose(FID);

s = size(v);

ts = v(1,:); % timestamp
t  = v(2:4,:); % translation [tx,ty,tz]

% Conversion from quaternions [qx qy qz qw] to rotation matrix and
% angle-axis
R = zeros(3,3,s(2));
u = zeros(3,s(2));
for i=1:s(2)
    r_xyz = v(5:7,i);
    % IMPORTANT:
    % Change sign of Y and Z rotations when correction is needed from Blender convention
    % In Blender the camera is pointing in the opposite direction to Z axis
    r_xyz(2:3) = -r_xyz(2:3);
    R(:,:,i) = RotationYPR( deg2rad( r_xyz(end:-1:1) ) ) ;
    u(:,i)   = angle_axis_rotation( R(:,:,i) );
end

gt = struct('ts', ts, 't',t, 'R',R, 'u',u);

end

