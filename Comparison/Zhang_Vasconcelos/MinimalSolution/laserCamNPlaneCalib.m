function T = laserCamNPlaneCalib(PIcam, L)

N_PLANES = size(PIcam,2);

if size(L,2) ~= N_PLANES
    error('The number of input lines and planes should be the equal');
end

plane_set = nchoosek(1:N_PLANES,3);

for i=1:size(plane_set,1)
    T{i} = laserCam3PlaneCalib(PIcam(:,plane_set(i,:)),L(:,plane_set(i,:)));
end