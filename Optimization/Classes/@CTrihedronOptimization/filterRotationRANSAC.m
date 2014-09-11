function obj = filterRotationRANSAC( obj )

% Mask of existing correspondences
mask_exist = obj.mask_LRF_V;
idxs_exist = find( mask_exist );

corresps = [ obj.cam_N ; obj.LRF_V ];
feedback = true;
[R0, inliers] = ransacFitTransNormals(corresps, obj.RANSAC_Rotation_threshold, feedback);
% R0 can be used as initial estimate for optimization process
obj.R0 = R0;

% Find indexes of outliers in complete set of 3xN correspondences
idxs_inliers  = idxs_exist( inliers );
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

% DEBUG: TODO: Refactor variables
% Complete set of correspondences
% all_mask = [co.thereis_line];
% all_idxs = find( all_mask );
% aux = [co.lab_line];
% all_label = zeros(1, length(all_mask));
% all_label(all_mask) = aux;
% all_inliers = all_idxs( inliers );
% fprintf('Selected %d inliers\n', length(inliers))
% all_outliers = setdiff( 1:length(all_mask), all_inliers );
% all_mask( all_outliers ) = false;
% all_label( all_outliers ) = 0;
% debug = 0;
% if debug
%     mask_lines = [co.thereis_line];
%     hF = figure('units','normalized','outerposition',[0 0 1 1]);
%     for ncorresp = 1:length(all_mask);
%         %             idx = 20;
%         %             ncorresp = all_inliers( idx );
%         nco = 1 + floor( ( ncorresp -1)/3 ); % Nr of CO observation
%         k = ncorresp - (nco-1)*3; % Nr of axis (X,Y,Z)
%         n_scan = check(nco).scan_frame;
%         n_img  = check(nco).cam_frame;
%         
%         % Plot
%         subplot(122)
%         cla
%         img_params = check(nco).img_params;
%         img = rgb2gray( imgs(n_img).I );
%         corner_calib_plot(img_params, [], img)
%         
%         subplot(121)
%         cla
%         hold on, title(sprintf('Scan segmentation in observation #%d',nco))
%         plot( scans(n_scan).xy(1,:), scans(n_scan).xy(2,:), '.k' ), axis equal
%         plotLIDARframe( );
%         col = 'rgb';
%         xyz = 'XYZ';
%         if ~mask_lines(ncorresp)
%             title( sprintf('NO OBSERVATION\nObservation #%d\nAxis %c',nco,xyz(k)) )
%         else
%             if all_mask(ncorresp)
%                 title( sprintf('INLIER\nObservation #%d\nAxis %c',nco,xyz(k)) )
%             else
%                 title( sprintf('OUTLIER\nObservation #%d\nAxis %c',nco,xyz(k)) )
%             end
%             pts = check(nco).scan_inPts{k};
%             lin = check(nco).scan_track(k).lin;
%             plot( pts(1,:), pts(2,:), [col(k),'.'] );
%             plotHomLineWin( lin, col(k) );
%         end
%         keyboard
%     end
end