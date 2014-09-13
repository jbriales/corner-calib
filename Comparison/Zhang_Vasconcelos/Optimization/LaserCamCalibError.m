%
% e = LaserCamCalibError(T,PIcam,L)
%
% Computes calibration error of a transformation T between a Camera and a 
% Laser reference frame, given a Plane in camera frame and a line in Laser 
% frame.
% Error is the euclidean distance between the plane and the line in the 
% dual space.
%
% INPUT
%   T     - 4x4 homogeneous transformation matrix from Camera to Laser 
%           frame
%   PIcam - 4xN matrix with columns as N plane coordinates in camera frame
%   L     - 6xN matrix with columns as N plucker line coordinates in laser 
%           frame
%

function e = LaserCamCalibError(T,PIcam,L)

N_PLANES = size(PIcam,2);

TP = [[T(1:3,1:3); -T(1:3,4).'*T(1:3,1:3)] [0;0;0;1]];

e = zeros(1,N_PLANES);
for n=1:N_PLANES
    PIlrs = TP*PIcam(:,n);
    PIlrs = PIlrs(1:3)/PIlrs(4);
    LD    = PluckerDual(L(:,n));
    e(n)  = DistancePointLine(PIlrs,LD);
end