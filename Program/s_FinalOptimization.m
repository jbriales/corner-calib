% Final optimization
co = co0;
% RANSAC selection of optimal correspondences
method = 2;
switch method
    case 1
        %% Filter correspondences checking error for each complete CO estimated rotation
        % Complete set of correspondences
        all_mask = [co.thereis_line];
        aux = [co.lab_line];
        all_label = zeros(1, length(all_mask));
        all_label(all_mask) = aux;
        all_n = [co.R_c_w];
        all_n = all_n(:,all_mask);
        all_l = cell2mat([co.l]);
        
        mask_completeCO = [co.complete];
        idxs_completeCO = find(mask_completeCO);
        fprintf('There are %d correspondences\n',sum(all_mask))
        fprintf('There are %d complete CO\n',sum(mask_completeCO))
        thres_rot = 3e-3;
        inliers  =  cell(1,sum(mask_completeCO));
        Ninliers = zeros(1,sum(mask_completeCO));
        for i=idxs_completeCO
            %         err = dot( all_n, gt.R_c_s(:,1:2)*all_l );
            err = dot( all_n, co(i).R_c_s(:,1:2)*all_l );
            %         hHist = figure; hist( err )
            med = median( abs( err ) );
            inliers{i}  = find( abs(err) < thres_rot );
            Ninliers(i) = length(inliers{i});
            fprintf('i=%d\tmedian=%e\t#inliers=%d\n',i,med,Ninliers(i))
            %         close( hHist )
        end
        [Ninliers_max,I] = max(Ninliers);
        fprintf('The selected rotation is the %d-th with %d inliers\n', I, Ninliers_max)
        inliers  = inliers{I};
        outliers = setdiff( 1:length(all_mask), inliers );
        all_mask( outliers ) = false;
        all_label( outliers ) = 0;
        
        % Set input data (correspondences)
        cont = 1;
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
                
                rot_input(cont).N = N;
                rot_input(cont).A_N = A_N;
                rot_input(cont).l = l;
                rot_input(cont).A_l = A_l;
                cont = cont + 1;
            end
        end
        
    
    case 2
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
        thres_ort = 100 * 7e-2;
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
        cont = 1;
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
                
                rot_input(cont).N = N;
                rot_input(cont).A_N = A_N;
                rot_input(cont).l = l;
                rot_input(cont).A_l = A_l;
                cont = cont + 1;
            end
        end
        
    case 3
        %% Filter complete CO with RANSAC and angular distances
        % RANSAC selection of optimal complete observations
        mask_complete = [co.complete];
        complete_idxs = find( mask_complete );
        all_rots = [co(mask_complete).R_c_s];
        all_rots = reshape( all_rots, 9, [] );
        [R, inliers] = ransacFilterRots(all_rots, 0.5, true);
        co = co(complete_idxs(inliers)); % TODO: Use the rest of information in some way
        fprintf('%d observations were kept after RANSAC on rotations\n',length(co));
        
        % Set input data (correspondences)
        for i=1:length(co)
            if co(i).complete
                N   = co(i).R_c_w;
                A_N = co(i).A_R_c_w;
                
                l   = cell2mat(co(i).l);
                A_l = diag(cell2mat(co(i).A_l));
            else
                occ = cellfun(@(x)~isempty(x), co(i).l);
                occ_= kron(ones(1,3),occ);
                occ_ = logical(occ_);
                
                N   = co(i).R_c_w(:,occ);
                A_N = co(i).A_R_c_w(occ_,occ_);
                
                l   = cell2mat(co(i).l(occ));
                A_l = diag(cell2mat(co(i).A_l(occ)));
            end
            rot_input(i).N = N;
            rot_input(i).A_N = A_N;
            rot_input(i).l = l;
            rot_input(i).A_l = A_l;
        end
end

% Average of calibration rotation matrices
if isfield(co,'R_c_s')
    R0 = sum( reshape([co.R_c_s],3,3,[]), 3 );
    [U,~,V] = svd( 0 );
    R0 = U*V';
else
    R0 = [ 0 -1  0
              0  0 -1
              1  0  0 ];
end

% s6_solveOptim
% s6_solveRotation % TODO: Change to optimRotation
[ R_c_s_w, cov_w, cov_eps_w, err_w, ~, ~ ] = optimRotation( rot_input, R0, true, all_label );
[ R_c_s_nw, cov_nw, cov_eps_nw, err_nw, ~, ~ ] = optimRotation( rot_input, R0, false, all_label );


% figure, hold on
% all_label( all_label==0 ) = [];
% rgb = 'rgb';
% subplot(221), hold on
% for i=1:length(all_label)
%     bar(i,err(i), rgb(all_label(i)))
% end
% subplot(223), hold on
% for i=1:length(all_label)
%     bar(i,W(i,i), rgb(all_label(i)))
% end
% subplot(222), hold on
% for i=1:length(all_label)
%     bar(i,err_(i), rgb(all_label(i)))
% end
% subplot(224), hold on
% for i=1:length(all_label)
%     bar(i,W_(i,i), rgb(all_label(i)))
% end

