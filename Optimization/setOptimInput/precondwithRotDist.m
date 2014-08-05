co = co0;
        
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

% Average of calibration rotation matrices
if isfield(co,'R_c_s')
    R_c_s = sum( reshape([co.R_c_s],3,3,[]), 3 );
    [U,~,V] = svd( R_c_s );
    R_c_s = U*V';
else
    R_c_s = [ 0 -1  0
              0  0 -1
              1  0  0 ];
end