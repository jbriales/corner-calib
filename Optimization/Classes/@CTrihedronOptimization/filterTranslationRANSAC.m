function obj = filterTranslationRANSAC( obj )

% Mask of existing correspondences
mask_exist = obj.mask_LRF_Q;
idxs_exist = find( mask_exist );

% Normalize P2 homogeneous lines to plane normals
N = obj.cam_L;
N = N ./ repmat( sqrt(sum(N.^2,1)) , 3,1);
corresps = [ N ; obj.LRF_Q ];
feedback = true;

s = 3;  % Minimum No of points needed to define rotation.
        
fittingfn = @computeTranslationLinear;
distfn    = @correspdist;
degenfn   = @isdegenerate;

[t0, inliers] = ransac(corresps, fittingfn, distfn, degenfn, s, obj.RANSAC_Translation_threshold, obj.debug_level);
   
% R0 can be used as initial estimate for optimization process
obj.t0 = t0;

% Find indexes of outliers in complete set of 3xN correspondences
idxs_inliers  = idxs_exist( inliers );
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

%% RANSAC auxiliary functions
function t = computeTranslationLinear( X )
L = X(1:3,:);
Q = X(4:5,:);

b = - dot( L, obj.R(:,1:2) * Q, 1 )';
A = L';
t = A \ b;
end

function [inliers, t] = correspdist(t, X, thres)
    
    all_n = X(1:3,:);
    all_q = X(4:5,:);
    Ncorresps = size(X,2);
    
    d = dot( all_n, obj.R(:,1:2)*all_q + repmat(t,1,Ncorresps) );
    
    inliers = find(abs(d) < thres);
end

function r = isdegenerate(X)
    % TODO: Degeneracy (rank) should be avoided to get good estimate
%     disp( X(1:3,:) )
%     fprintf('Det  = %f\n', abs( det(X(1:3,:)) ));
%     fprintf('Cond = %f\n', cond( X(1:3,:) ));
%     r = false;
    r = ( abs(det(X(1:3,:))) < 1e-3 ); % Best value seems to be 1e-2
end

end
