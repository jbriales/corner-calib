function obj = filterRotation_R0( obj, R0 )

% Mask of existing correspondences
mask_exist = obj.mask_LRF_V;
idxs_exist = find( mask_exist );

% Compute error and check inliers
res = obj.FErr_Orthogonality( R0 );
inliers = find( abs(res) < obj.RANSAC_Rotation_threshold );

obj.R0 = R0;

% Find indexes of outliers in complete set of 3xN correspondences
idxs_inliers  = idxs_exist( inliers ); %#ok<FNDSB>
idxs_outliers = setdiff( idxs_exist, idxs_inliers );

%% Assign outlier tag to each observation obs(i)
Nobs = length( obj.obs );
map = reshape( 1:3*Nobs, 3,Nobs );

% Assign 1 to elements which are outliers
mask_RANSAC_R_outliers = false(1,length(mask_exist));
mask_RANSAC_R_outliers( idxs_outliers ) = true;
for i=1:Nobs
    obj.obs(i).is_R_outlier = mask_RANSAC_R_outliers( map(:,i) );
end

end