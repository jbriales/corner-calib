function [best_T, best_inliers, best_sample] = laserCamLinearRobust(NPI,LrsPts,threshold)

for i=1:length(LrsPts)
    if size(LrsPts{i},1) == 2
        LrsPts{i} = [LrsPts{i}; ones(1,size(LrsPts{i},2))];
    end
end 

N_PLANES = size(NPI,2);

PIcam = [-NPI./((ones(3,1)*sqrt(sum(NPI.^2))).^2); ones(1,N_PLANES)];

plane_set = nchoosek(1:N_PLANES,5);

min_cost = inf;
error = nan(1,5);
for i = 1:size(plane_set,1)
    T = ZhangAlgorithm(NPI(:,plane_set(i,:)),LrsPts(plane_set(i,:)));
    for j=1:N_PLANES
        error(j) = mean(laserCamRes(LrsPts{j},PIcam(:,j),T));
    end
    inliers = union(find(error < threshold),plane_set(i,:));
    cost    = sum(error(inliers).^2) + threshold^2*(size(error,2)-length(inliers));
    if cost < min_cost
        min_cost     = cost;
        best_inliers = inliers;
        best_T       = T;
        best_sample  = plane_set(i,:);
    end

end