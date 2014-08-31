co = co0;

%% Filter correspondences with RANSAC in all correspondences
% Complete set of correspondences
all_mask = [co.thereis_line];
all_idxs = find( all_mask );
aux = [co.lab_line];
all_label = zeros(1, length(all_mask));
all_label(all_mask) = aux;
all_n = [co.R_c_w];
all_n = all_n(:,all_mask);
all_l = cell2mat([co.l]);

corresps = [ all_n ; all_l ];
%         thres_ort = 1e-3;
thres_ort = 1e-1;
feedback = true;
[R, inliers] = ransacFitTransNormals(corresps, thres_ort, feedback);
all_inliers = all_idxs( inliers );
fprintf('Selected %d inliers\n', length(inliers))
all_outliers = setdiff( 1:length(all_mask), all_inliers );
all_mask( all_outliers ) = false;
all_label( all_outliers ) = 0;

debug = 0;
if debug
    mask_lines = [co.thereis_line];
    hF = figure('units','normalized','outerposition',[0 0 1 1]);
    for ncorresp = 1:length(all_mask);
        %             idx = 20;
        %             ncorresp = all_inliers( idx );
        nco = 1 + floor( ( ncorresp -1)/3 ); % Nr of CO observation
        k = ncorresp - (nco-1)*3; % Nr of axis (X,Y,Z)
        n_scan = check(nco).scan_frame;
        n_img  = check(nco).cam_frame;
        
        % Plot
        subplot(122)
        cla
        img_params = check(nco).img_params;
        img = rgb2gray( imgs(n_img).I );
        corner_calib_plot(img_params, [], img)
        
        subplot(121)
        cla
        hold on, title(sprintf('Scan segmentation in observation #%d',nco))
        plot( scans(n_scan).xy(1,:), scans(n_scan).xy(2,:), '.k' ), axis equal
        plotLIDARframe( );
        col = 'rgb';
        xyz = 'XYZ';
        if ~mask_lines(ncorresp)
            title( sprintf('NO OBSERVATION\nObservation #%d\nAxis %c',nco,xyz(k)) )
        else
            if all_mask(ncorresp)
                title( sprintf('INLIER\nObservation #%d\nAxis %c',nco,xyz(k)) )
            else
                title( sprintf('OUTLIER\nObservation #%d\nAxis %c',nco,xyz(k)) )
            end
            pts = check(nco).scan_inPts{k};
            lin = check(nco).scan_track(k).lin;
            plot( pts(1,:), pts(2,:), [col(k),'.'] );
            plotHomLineWin( lin, col(k) );
        end
        keyboard
    end
end

% Set input data (correspondences)
count = 1;
for i=1:numel(all_mask)/3
    k = 1 + (i-1)*3;
    
    occ = all_mask(k:k+2);
    
    if any(occ)
        occ_= kron(ones(1,3),occ);
        occ_ = logical(occ_);
        
        N   = co(i).R_c_w(:,occ);
        A_N = co(i).A_R_c_w(occ_,occ_);
        
        l   = cell2mat(co(i).l(occ));
        A_l = diag(cell2mat(co(i).A_l(occ)));
        
        rot_input(count).N = N;
        rot_input(count).A_N = A_N;
        rot_input(count).l = l;
        rot_input(count).A_l = A_l;
        
        % For Monte Carlo simulation
        J_eps_R = MonteCarlo.manDiffRotLog( co(i).R_c_w ); 
        rot_input(count).A_eps = J_eps_R * co(i).A_R_c_w * J_eps_R';
        
        % Increase counter
        count = count + 1;
    end
end

% Average of calibration rotation matrices
if isfield(co,'R_c_s')
    R0 = sum( reshape([co.R_c_s],3,3,[]), 3 );
    [U,~,V] = svd( R0 );
    R0 = U*V';
else
    R0 = [ 0 -1  0
              0  0 -1
              1  0  0 ];
end