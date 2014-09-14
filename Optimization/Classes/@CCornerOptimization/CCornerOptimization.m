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
        cam_L
        
        LRF_Q
        
        Rt0
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
        % For Rotation and translation together
        function residual = FErr_2D_LineDistance( obj, Rt )
            % Compute error vector for observations data and certain Rt
            L = obj.cam_L;
            Q = obj.LRF_Q;
            % Check normalized
            R = Rt(1:3,1:3);
            t = Rt(1:3,4);
            Ncorr = size( Q,2 );
            P = hnormalise( obj.K * ( R(:,1:2) * Q + repmat( t,1,Ncorr ) ) );
            residual = dot( L, P, 1 )';
        end
        function jacobian = FJac_2D_LineDistance( obj, Rt )
            % Compute jacobian of error vector wrt R for observations data and certain R
            L = obj.cam_L;
            Q = obj.LRF_Q;
            
            R = Rt(1:3,1:3);
            t = Rt(1:3,4);
                  
            Ncorr = size( Q,2 );
            jacobian = zeros(Ncorr,6);
            for i=1:Ncorr
                % Homogeneous 3D points in camera frame (pixel units)
                %             Ph = obj.K * ( R(:,1:2) * Q + repmat( t,1,Ncorr ) );
                
                l = L(:,i);
                q = Q(:,i);
                ph = obj.K * ( R(:,1:2) * q + t );
                J_hnormalise = 1/ph(3)^2 * [ ph(3) 0 -ph(1)
                                             0 ph(3) -ph(2)
                                               zeros(1,3)   ];
                J_ph_Rt = obj.K * [ -skew(R(:,1:2)*q) , eye(3) ];
                
                jacobian(i,:) = l' * J_hnormalise * J_ph_Rt;
            end
        end
        function weights  = FWeights_2D_LineDistance( obj, Rt )
            residual = obj.FErr_2D_LineDistance( Rt );
            residual = reshape( residual, 3, [] );
            N = size(residual,2);
            weights  = dot( residual, residual, 2 ) / N;
            weights  = kron( eye(N), diag(1./weights) );
        end
        
        function [R,t] = optimizeRt_NonWeighted( obj )
            Fun = @(Rt) deal( obj.FErr_2D_LineDistance( Rt ) , obj.FJac_2D_LineDistance( Rt ) );
            Rt = obj.optimize( Fun, obj.Rt0, 'SE(3)', false );
            R = Rt(1:3,1:3);
            t = Rt(1:3,4);
        end
        function [R,t] = optimizeRt_Weighted( obj )
            Fun = @(Rt) deal( obj.FErr_2D_LineDistance( Rt ) ,...
                              obj.FJac_2D_LineDistance( Rt ) ,...
                              obj.FWeights_2D_LineDistance( Rt ) );
            Rt = obj.optimize( Fun, obj.Rt0, 'SE(3)', true );
            R = Rt(1:3,1:3);
            t = Rt(1:3,4);
        end
        function [R,t] = optimizeRt_ConstWeighted( obj )
            [R,t] = obj.optimizeRt_NonWeighted;
            Rt0 = [R,t];
            weights = obj.FWeights_2D_LineDistance( Rt0 );
            
            Fun = @(Rt) deal( obj.FErr_2D_LineDistance( Rt ) ,...
                              obj.FJac_2D_LineDistance( Rt ) ,...
                              weights );
            Rt = obj.optimize( Fun, Rt0, 'SE(3)', true );
            R = Rt(1:3,1:3);
            t = Rt(1:3,4);
        end
        function [R,t] = optimizeRt_PreWeighted( obj )
            [R,t] = obj.optimizeRt_NonWeighted;
            Rt0 = [R,t];
            
            Fun = @(Rt) deal( obj.FErr_2D_LineDistance( Rt ) ,...
                              obj.FJac_2D_LineDistance( Rt ) ,...
                              obj.FWeights_2D_LineDistance( Rt ) );
            Rt = obj.optimize( Fun, Rt0, 'SE(3)', true );
            R = Rt(1:3,1:3);
            t = Rt(1:3,4);
        end
        
        function h = plotRotationCostFunction( obj, R, t )
            gv  = obj.get_plot_gv( obj.plot_dist_R );
            
            FE = @(R)obj.FErr_2D_LineDistance( [R,t] );
            W  = obj.FWeights_2D_LineDistance( [R,t] );
            Fx = @(R,inc) expmap( inc ) * R;
            x0 = R;
            h  = obj.plotCostFunction( gv, W, FE, Fx, x0 );
        end
        
        function h = plotTranslationCostFunction( obj, R, t )
            gv  = obj.get_plot_gv( obj.plot_dist_t );
            
            FE = @(t)obj.FErr_2D_LineDistance( [R,t] );
            W  = obj.FWeights_2D_LineDistance( [R,t] );
            Fx = @(t,inc) t + inc;
            x0 = t;
            h  = obj.plotCostFunction( gv, W, FE, Fx, x0 );
        end
        
        plot( obj )
        
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
        % For translation
        function cam_L = get.cam_L( obj )
            cam_L = [obj.obs.cam_l];
            % Mask with existing measures and not-inliers
%             mask_valid = obj.mask_LRF_Q & (~obj.mask_RANSAC_t_outliers);
%             cam_L = cam_L(:, mask_valid);
            cam_L = cell2mat( cam_L );
        end
%         
%         function cam_reprN = get.cam_reprN( obj )
%             cam_reprN = [obj.obs.cam_reprN];
%             % Mask with existing measures and not-inliers
%             mask_valid = obj.mask_LRF_Q & (~obj.mask_RANSAC_t_outliers);
%             cam_reprN = cam_reprN(:, mask_valid);
%         end
%         
        function LRF_Q = get.LRF_Q( obj )
            LRF_Q = [obj.obs.LRF_q];
            % Remove non-observed data: Not necessary really (already empty)
            % Remove outliers: Set as empty observations
%             mask_valid = obj.mask_LRF_Q & (~obj.mask_RANSAC_t_outliers);
%             LRF_Q( ~mask_valid ) = [];
            % Convert from cell with empty elements to dense matrix
            LRF_Q = cell2mat( LRF_Q );
        end
        
%         function mask_LRF_V = get.mask_LRF_V( obj )
%             mask_LRF_V = [obj.obs.thereis_LRF_v];
%         end
%         function mask_LRF_Q = get.mask_LRF_Q( obj )
%             mask_LRF_Q = [obj.obs.thereis_LRF_q];
%         end

        function Rt0 = get.Rt0( obj )
            Rt0 = [ obj.R0, obj.t0 ];
        end

    end
    
end

