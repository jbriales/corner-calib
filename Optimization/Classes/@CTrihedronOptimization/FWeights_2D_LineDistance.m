function weights = FWeights_2D_LineDistance( obj, R, t )
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
    mask_L = obj.obs(i).thereis_LRF_q & (~obj.obs(i).is_t_outlier);
    
    % Homogeneous lines in calibrated (m) space
    L = obj.obs(i).cam_l( :, mask_L );
    % Convert L to uncalibrated (pixel) space and normalize line:
    L = obj.K' \ L;
    L = L ./ repmat( sqrt(sum(L(1:2,:).^2,1)), 3,1 );
    
    % Correlated covariance of calibrated homogeneous lines
    A_L = obj.obs(i).cam_A_l;
    % Convert covariance to uncalibrated (pixel) space
    Jac_pix_m = kron(eye(3), inv(obj.K'));
    A_L = Jac_pix_m * A_L * Jac_pix_m';
    A_L = mat2cell(A_L,[3 3 3],[3 3 3]);
    A_L = A_L(mask_L,mask_L);
    A_L = cell2mat(A_L);
    
    mask_Q = obj.obs(i).thereis_LRF_q & (~obj.obs(i).is_t_outlier);
    Q = obj.obs(i).LRF_q;
    Q( ~mask_Q ) = [];
    Q = cell2mat( Q );
    
    % Could be used non-minimal covariance
    % Input LRF_A_v is a 1x3 cell array (diagonal elements)
    A_Q = obj.obs(i).LRF_A_q;
    % Remove outlier of existing part of A_V
    A_Q(~mask_Q) = [];
    A_Q = blkdiag( A_Q{:} );
    
    % TODO: Modify
    obsSize  = size( L, 2 );
    if obsSize > 0
        % Compute all jacobians wrt L_XYZ
        Ph = hnormalise( obj.K * ( R(:,1:2) * Q + repmat(t,1,obsSize) ) );
        Ph = mat2cell( Ph, 3, ones(1,obsSize) );
        JN = blkdiag( Ph{:} )'; % Derivative of e_XYZ wrt L_XYZ
        
        % Compute all jacobians wrt q_XYZ
        JQ = cell(1,obsSize);
        for k=1:obsSize
            ph = Ph{k};
            J_hnormalise = 1/ph(3)^2 * [ ph(3) 0 -ph(1)
                                         0 ph(3) -ph(2)
                                            zeros(1,3)   ];
            JQ{k} = L(:,k)' * J_hnormalise * obj.K * R(:,1:2);
        end
        JQ = blkdiag( JQ{:} );

        A_i = JN * A_L * JN' + JQ * A_Q * JQ';
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