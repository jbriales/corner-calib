function obj = filterTranslation_t0( obj, R, t0 )

% Mask of existing correspondences
mask_exist = obj.mask_LRF_Q;
idxs_exist = find( mask_exist );

res = obj.FErr_3D_PlaneDistance( R, t0 );
inliers = find( abs(res) < obj.RANSAC_Translation_threshold );

% t0 can be used as initial estimate for optimization process
obj.t0 = t0;

% Find indexes of outliers in complete set of 3xN correspondences
idxs_inliers  = idxs_exist( inliers ); %#ok<FNDSB>
idxs_outliers = setdiff( idxs_exist, idxs_inliers );

%% Assign outlier tag to each observation obs(i)
Nobs = length( obj.obs );
map = reshape( 1:3*Nobs, 3,Nobs );

% Assign 1 to elements which are outliers
mask_RANSAC_t_outliers = false(1,length(mask_exist));
mask_RANSAC_t_outliers( idxs_outliers ) = true;
for i=1:Nobs
    obj.obs(i).is_t_outlier = mask_RANSAC_t_outliers( map(:,i) );
end

end
