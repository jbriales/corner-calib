classdef CTrihedronOptimization < handle
    %CTrihedronOptimization Class to store observations and later process
    %and optimize them
    %   Detailed explanation goes here
    
    properties
        obs     % Array of CTrihedronObservation
        
        % Auxiliar variables
        K   % Camera intrinsic matrix (for translation optimization)
        
        RANSAC_Rotation_threshold % Threshold for rotation error function
        RANSAC_Translation_threshold % Threshold for translation error function
        debug_level % Verbose level when optimizing
        
        maxIters    % Max number of iterations in LM optimization
    end
    
    properties (SetAccess=private)
        % Optimized variables
        R0
        % R % R is not a property to avoid confusion with results of
        % different methods
        t0
        % t % t is not a property to avoid confusion with results of
        % different methods
    end
    
    properties (Dependent) % Array values collected from set of observations
        cam_N
        cam_L
        cam_reprN
        
        LRF_V
        LRF_L
        LRF_Q
        
        mask_LRF_V
        mask_LRF_Q
        mask_RANSAC_R_outliers
        mask_RANSAC_t_outliers
    end
    
    methods
        %% Constructor
        function obj = CTrihedronOptimization( K, RANSAC_Rotation_threshold, RANSAC_Translation_threshold, debug_level, maxIters )
            obj.obs = CTrihedronObservation.empty(1,0);
            
            if ~exist('RANSAC_Rotation_threshold','var')
                RANSAC_Rotation_threshold = 1e-1;
            end
            obj.RANSAC_Rotation_threshold = RANSAC_Rotation_threshold;
            
            if ~exist('RANSAC_Translation_threshold','var')
                RANSAC_Translation_threshold = 1e-3;
            end
            obj.RANSAC_Translation_threshold = RANSAC_Translation_threshold;
            
            if ~exist('debug_level','var')
                debug_level = 2;
            end
            obj.debug_level = debug_level;
            
            if ~exist('maxIters','var')
                maxIters = 50;
            end
            obj.maxIters = maxIters;
            
            obj.K = K;
            
            % Commented because is dependent now
%             obj.mask_RANSAC_R_outliers = []; % Initialize as empty
        end
        
        %% Load new observation
        function obj = stackObservation( obj, obs )
            if isa( obs, 'CTrihedronObservation' )
                obj.obs(end+1) = obs;
            else
                error('The loaded observation is not of class ''CTrihedronObservation''');
            end
        end
        
        %% Filter correspondences with RANSAC
        obj = filterRotationRANSAC( obj )
        obj = filterTranslationRANSAC( obj, R )
        % Functions to set all outlier masks to false
        function obj = resetRotationRANSAC( obj )
            Nobs = length( obj.obs );
            for i=1:Nobs
                obj.obs(i).is_R_outlier = false(1,3);
            end
        end
        function obj = resetTranslationRANSAC( obj )
            Nobs = length( obj.obs );
            for i=1:Nobs
                obj.obs(i).is_t_outlier = false(1,3);
            end
        end
        
        %% Compute initial estimate
        function obj = setInitialRotation( obj, R0 )
            obj.R0 = R0;
        end
        
        function obj = setInitialTranslation( obj, t0 )
            obj.t0 = t0;
        end
        
        function obj = computeTranslationLinear( obj )
            L = obj.cam_L;
            q = obj.LRF_Q;
            
            b = - dot( L, obj.R(:,1:2) * q, 1 )';
            A = L';
            t_lin = A \ b;
            
            obj.t0 = t_lin;
        end
        % TODO: Automate from complete poses
        
        %% Optimization functions and methods
        % For Rotation
        function residual = FErr_Orthogonality( obj, R )
            % Compute error vector for observations data and certain R
            n_cam = obj.cam_N;
            v_LRF = obj.LRF_V;
            residual = dot( n_cam, R(1:3,1:2) * v_LRF, 1 )';
        end
        function jacobian = FJac_Orthogonality( obj, R )
            % Compute jacobian of error vector wrt R for observations data and certain R
            n_cam = obj.cam_N;
            v_LRF = obj.LRF_V;
            jacobian = cross( R(1:3,1:2) * v_LRF, n_cam, 1 )';
        end
        weights = FWeights_Orthogonality( obj, R )
        
        R = optimizeRotation_NonWeighted( obj )
        R = optimizeRotation_Weighted( obj )
        R = optimizeRotation_DiagWeighted( obj )
        J = FJac_optimizeRotation_NonWeighted( obj, R )
        
        % TODO: Implement optimizeRotation_Covariance from optimRotation.m
        
        h = plotRotationCostFunction( obj, R )
        
        plot( obj )
        
        % For Translation
        function residual = FErr_3D_PlaneDistance( obj, R, t )
            % Compute error vector for observations data, R and t
            % N - (3x...) 3D normals to reprojection planes from camera center
            % through image lines
            % Q - (2x...) 2D LRF intersection points
            L = obj.cam_L;
            N = L ./ repmat( sqrt(sum(L.^2,1)), 3,1 ); % Normalize vectors as plane normals
            Q = obj.LRF_Q;
            
            s = size(N,2); % Number of correspondences
            
            residual = dot( N, R(:,1:2)*Q + repmat(t,1,s), 1)';
        end
        function jacobian = FJac_3D_PlaneDistance( obj, R, t )
            % Compute jacobian of error vector wrt R for observations data and certain R
            L = obj.cam_L;
            N = L ./ repmat( sqrt(sum(L.^2,1)), 3,1 ); % Normalize vectors as plane normals
            
            jacobian = N';
        end
        weights = FWeights_3D_PlaneDistance( obj, R, t )
        
        t = optimizeTranslation_2D_NonWeighted( obj, R )
        t = optimizeTranslation_2D_Weighted( obj, R )
        t = optimizeTranslation_3D_NonWeighted( obj, R )
        t = optimizeTranslation_3D_Weighted( obj, R )
        
        h = plotTranslation_3D_CostFunction( obj, R, t )
        
        
        %% Get-functions
        % For rotation
        function cam_N = get.cam_N( obj )
            cam_N = [obj.obs.cam_R_c_w];
            % Mask with existing measures and not-inliers
            mask_valid = obj.mask_LRF_V & (~obj.mask_RANSAC_R_outliers);
            cam_N = cam_N(:, mask_valid);
        end
        
        function LRF_V = get.LRF_V( obj )
            LRF_V = [obj.obs.LRF_v];
            % Remove non-observed data: Not necessary really (already empty)
            % Remove outliers: Set as empty observations
            mask_valid = obj.mask_LRF_V & (~obj.mask_RANSAC_R_outliers);
            LRF_V( ~mask_valid ) = [];
            % Convert from cell with empty elements to dense matrix
            LRF_V = cell2mat( LRF_V );
        end
        
        % For translation
        function cam_L = get.cam_L( obj )
            cam_L = [obj.obs.cam_l];
            % Mask with existing measures and not-inliers
            mask_valid = obj.mask_LRF_Q & (~obj.mask_RANSAC_t_outliers);
            cam_L = cam_L(:, mask_valid);
        end
        
        function cam_reprN = get.cam_reprN( obj )
            cam_reprN = [obj.obs.cam_reprN];
            % Mask with existing measures and not-inliers
            mask_valid = obj.mask_LRF_Q & (~obj.mask_RANSAC_t_outliers);
            cam_reprN = cam_reprN(:, mask_valid);
        end
        
        function LRF_Q = get.LRF_Q( obj )
            LRF_Q = [obj.obs.LRF_q];
            % Remove non-observed data: Not necessary really (already empty)
            % Remove outliers: Set as empty observations
            mask_valid = obj.mask_LRF_Q & (~obj.mask_RANSAC_t_outliers);
            LRF_Q( ~mask_valid ) = [];
            % Convert from cell with empty elements to dense matrix
            LRF_Q = cell2mat( LRF_Q );
        end
        
        function mask_LRF_V = get.mask_LRF_V( obj )
            mask_LRF_V = [obj.obs.thereis_LRF_v];
        end
        function mask_LRF_Q = get.mask_LRF_Q( obj )
            mask_LRF_Q = [obj.obs.thereis_LRF_q];
        end
        function mask_RANSAC_R_outliers = get.mask_RANSAC_R_outliers( obj )
            mask_RANSAC_R_outliers = [obj.obs.is_R_outlier];
        end
        function mask_RANSAC_t_outliers = get.mask_RANSAC_t_outliers( obj )
            mask_RANSAC_t_outliers = [obj.obs.is_t_outlier];
        end
        
        % More functions
        function disp_N_R_inliers( obj )
            mask_valid = obj.mask_LRF_V & (~obj.mask_RANSAC_R_outliers);
            N = sum( mask_valid );
            Ncorresp = sum( obj.mask_LRF_V );
            fprintf('The number of R inliers is %d of %d\n',N,Ncorresp);
        end
        function disp_N_t_inliers( obj )
            mask_valid = obj.mask_LRF_Q & (~obj.mask_RANSAC_t_outliers);
            N = sum( mask_valid );
            Ncorresp = sum( obj.mask_LRF_Q );
            fprintf('The number of t inliers is %d of %d\n',N,Ncorresp);
        end
    end
    
end

