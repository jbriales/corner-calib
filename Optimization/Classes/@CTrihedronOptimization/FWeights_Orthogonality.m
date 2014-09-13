function weights = FWeights_Orthogonality( obj, R )
% Compute covariance matrix of error vector for observations
% data and certain R and invert it for using as weight matrix W
N_obs = length(obj.obs);
W = cell(1, N_obs);

ort = [ 0 -1
    1  0 ];

R12 = R(1:3,1:2); % Only 2 columns are used

for i=1:N_obs
    % Take observation elements which exist and are not
    % outliers
    mask_N = obj.obs(i).thereis_LRF_v & (~obj.obs(i).is_R_outlier);
    
    % Normals to planes
    N = obj.obs(i).cam_R_c_w( :, mask_N );
    
    % Correlated covariance of normals to planes
    A_N = obj.obs(i).cam_A_R_c_w;
    A_N = mat2cell(A_N,[3 3 3],[3 3 3]);
    A_N = A_N(mask_N,mask_N);
    A_N = cell2mat(A_N);
    
    mask_V = obj.obs(i).thereis_LRF_v & (~obj.obs(i).is_R_outlier);
    V = obj.obs(i).LRF_v;
    V( ~mask_V ) = [];
    V = cell2mat( V );
    
    % Could be used non-minimal covariance
    % Input LRF_A_v is a 1x3 cell array (diagonal elements)
    A_V = obj.obs(i).LRF_A_v;
    % Remove outlier of existing part of A:V
    A_V(~mask_V) = [];
    A_V = diag(cell2mat(A_V));
    
    obsSize  = size( N, 2 );
    RL = mat2cell( R12*V, 3, ones(1,obsSize) );
    JN = blkdiag( RL{:} )'; % Derivative of e_XYZ wrt N_XYZ
    Ja = diag( dot( N, R12*ort*V ) );
    A_i = JN * A_N * JN' + Ja * A_V * Ja';
    if any(A_i(:))
        W{i} = pinv( A_i );
    else
        W{i} = eye(length(A_i))
    end
end
weights = blkdiag( W{:} );
end