classdef CCornerOptimization < handle & CBaseOptimization
    %CCornerOptimization Class to store observations and later process
    %and optimize them
    %   Detailed explanation goes here
    
    properties
        obs     % Array of CCornerObservation
        
        % Auxiliar variables
        K   % Camera intrinsic matrix (for translation optimization)
        
    end
    
    properties (SetAccess=private)
        % Optimized variables
        R0
        t0
    end
    
    properties (Dependent) % Array values collected from set of observations

    end
    
    methods
        %% Constructor
        function obj = CCornerOptimization( K, ...
                debug_level,maxIters, minParamChange, minErrorChange ) % COptimization inputs
            % Set optimization parameters through parent class
            obj = obj@CBaseOptimization( debug_level, maxIters, minParamChange, minErrorChange );

            obj.obs = CCornerObservation.empty(1,0);
            
            obj.K = K;
        end
        
        %% Load new observation
        function obj = stackObservation( obj, obs )
            obj.obs(end+1) = obs;
        end
                
        %% Compute initial estimate
        function obj = setInitialRotation( obj, R0 )
            obj.R0 = R0;
        end
        
        function obj = setInitialTranslation( obj, t0 )
            obj.t0 = t0;
        end
        
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
%         % For rotation
%         function cam_N = get.cam_N( obj )
%             cam_N = [obj.obs.cam_R_c_w];
%             % Mask with existing measures and not-inliers
%             mask_valid = obj.mask_LRF_V & (~obj.mask_RANSAC_R_outliers);
%             cam_N = cam_N(:, mask_valid);
%         end
%         
%         function LRF_V = get.LRF_V( obj )
%             LRF_V = [obj.obs.LRF_v];
%             % Remove non-observed data: Not necessary really (already empty)
%             % Remove outliers: Set as empty observations
%             mask_valid = obj.mask_LRF_V & (~obj.mask_RANSAC_R_outliers);
%             LRF_V( ~mask_valid ) = [];
%             % Convert from cell with empty elements to dense matrix
%             LRF_V = cell2mat( LRF_V );
%         end
%         
%         % For translation
%         function cam_L = get.cam_L( obj )
%             cam_L = [obj.obs.cam_l];
%             % Mask with existing measures and not-inliers
%             mask_valid = obj.mask_LRF_Q & (~obj.mask_RANSAC_t_outliers);
%             cam_L = cam_L(:, mask_valid);
%         end
%         
%         function cam_reprN = get.cam_reprN( obj )
%             cam_reprN = [obj.obs.cam_reprN];
%             % Mask with existing measures and not-inliers
%             mask_valid = obj.mask_LRF_Q & (~obj.mask_RANSAC_t_outliers);
%             cam_reprN = cam_reprN(:, mask_valid);
%         end
%         
%         function LRF_Q = get.LRF_Q( obj )
%             LRF_Q = [obj.obs.LRF_q];
%             % Remove non-observed data: Not necessary really (already empty)
%             % Remove outliers: Set as empty observations
%             mask_valid = obj.mask_LRF_Q & (~obj.mask_RANSAC_t_outliers);
%             LRF_Q( ~mask_valid ) = [];
%             % Convert from cell with empty elements to dense matrix
%             LRF_Q = cell2mat( LRF_Q );
%         end
        
%         function mask_LRF_V = get.mask_LRF_V( obj )
%             mask_LRF_V = [obj.obs.thereis_LRF_v];
%         end
%         function mask_LRF_Q = get.mask_LRF_Q( obj )
%             mask_LRF_Q = [obj.obs.thereis_LRF_q];
%         end

    end
    
end

