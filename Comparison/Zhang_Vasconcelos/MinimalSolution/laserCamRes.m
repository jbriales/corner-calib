% residue of calibration Tcam2lrs for given lrs-camera measurements 
% (Plrs,Tp2cam) lrs-camera
% 
%
% INPUTS:
%   Plrs     - 2xN matrix with 2D coordinates of N lrs point measurements
%   PIcam    - 4x1 matrix with projective coordinates of N planes in camera 
%              reference frame 
%   Tcam2lrs - 4x4 transformation matrix between camera and lrs reference
%              frames (lrs-camera calibration)


function r = laserCamRes(Plrs,PIcam,Tcam2lrs)

% plane in lrs reference frame
PIlrs  = Tcam2lrs.' \ PIcam;

% lrs measurement lines
L(1:3,:) = [Plrs(1:2,:); ones(1,size(Plrs,2))];
for j=1:size(Plrs,2)
    L(4:6,j) = cross(Plrs(:,j),[0;0;1]);
end

% intersect lines and planes
Pest = IntersectionLinePlane(PIlrs, L);

% get residue
r = sqrt(sum((Plrs(1:2,:) - Pest(1:2,:)).^2));
