function [best_T, best_inliers, best_sample] = laserCamMSACCalib(PIcam,L,threshold)

N_PLANES = size(PIcam,2);

plane_set = nchoosek(1:N_PLANES,3);

min_cost = inf;
for i = 1:size(plane_set,1)
    T_set = laserCam3PlaneCalib(PIcam(:,plane_set(i,:)), L(:,plane_set(i,:)));
    
    error = nan(1,N_PLANES);
    for m = 1:size(T_set,3)
        error(m,:) = LaserCamCalibError(T_set(:,:,m), PIcam, L);
        inliers    = union(find(error(m,:) < threshold),plane_set(i,:));
        cost       = sum(error(m,inliers).^2) + threshold^2*(size(error,2)-length(inliers));
        if cost < min_cost
            min_cost     = cost;
            best_inliers = inliers;
            best_T       = T_set(:,:,m);
            best_sample  = plane_set(i,:);
        end
    end
end