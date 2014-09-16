function weights = FWeights_3D_PlaneDistance( obj, R, t )
% Compute covariance matrix of error vector for observations
% data and certain R and t and invert it for using as weight matrix W
N_obs = obj.Nobs;
W = cell(1, N_obs);

ort = [ 0 -1
    1  0 ];

R12 = R(1:3,1:2); % Only 2 columns are used

for i=1:N_obs
    % Take observation elements which exist and are not
    % outliers
    mask_N = obj.obs(i).thereis_LRF_q & (~obj.obs(i).is_t_outlier);
    
    % Normals to reprojection planes
    N = obj.obs(i).cam_reprN( :, mask_N );
    
    % Correlated covariance of normals to planes
    A_N = obj.obs(i).cam_A_reprN;
    A_N = mat2cell(A_N,[3 3 3],[3 3 3]);
    A_N = A_N(mask_N,mask_N);
    A_N = cell2mat(A_N);
    
    mask_Q = obj.obs(i).thereis_LRF_q & (~obj.obs(i).is_t_outlier);
    Q = obj.obs(i).LRF_q;
    Q( ~mask_Q ) = [];
    Q = cell2mat( Q );
    
    % Could be used non-minimal covariance
    % Input LRF_A_v is a 1x3 cell array (diagonal elements)
    A_Q = obj.obs(i).LRF_A_q;
    % Remove outlier of existing part of A:V
    A_Q(~mask_Q) = [];
    A_Q = blkdiag( A_Q{:} );
    
    % TODO: Modify
    obsSize  = size( N, 2 );
    if obsSize > 0
        % Compute all jacobians wrt N_XYZ
        Rq_t = R12*Q + repmat(t,1,obsSize );
        Rq_t = mat2cell( Rq_t, 3, ones(1,obsSize) );
        JN = blkdiag( Rq_t{:} )'; % Derivative of e_XYZ wrt N_XYZ
        % Compute all jacobians wrt q_XYZ
        nR = N' * R12;
        nR = mat2cell( nR', 2, ones(1,obsSize) );
        JQ = blkdiag( nR{:} )';
        A_i = JN * A_N * JN' + JQ * A_Q * JQ';
        if any(A_i(:))
            W{i} = pinv( A_i );
        else
            W{i} = eye(length(A_i));
        end
    else
        W{i} = []; % Observation does not appear
    end
end
weights = blkdiag( W{:} );
end