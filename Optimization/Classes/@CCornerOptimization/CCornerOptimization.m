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
        % Arrays with central and virtual correspondences
        cam_L
        LRF_Q
        
        % Get arrays with central correspondences only
        cam_C_L
        LRF_C_Q
        
        Rt0
        
        Ncorresp
        Ncorresp_C
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
            if ~isempty( obs )
                obj.obs(end+1) = obs;
            end
        end
                
        %% Compute initial estimate
        function obj = setInitialRotation( obj, R0 )
            obj.R0 = R0;
        end
        
        function obj = setInitialTranslation( obj, t0 )
            obj.t0 = t0;
        end
        
        %% Optimization functions and methods
        % For coupled rotation and translation
        % Functions with central and virtual correspondences
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
        function H = FHes_2D_LineDistance( obj, Rt )
            % Linear approximation to hessian in LM method
            jacobian = obj.FJac_2D_LineDistance( Rt );
            weights  = obj.FWeights_2D_LineDistance( Rt );
            H = jacobian' * weights * jacobian;
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
        
        % Functions with central correspondences only
        function residual = FErr_C_2D_LineDistance( obj, Rt )
            % Compute error vector for observations data and certain Rt
            L = obj.cam_C_L;
            Q = obj.LRF_C_Q;
            % Check normalized
            R = Rt(1:3,1:3);
            t = Rt(1:3,4);
            Ncorr = size( Q,2 );
            P = hnormalise( obj.K * ( R(:,1:2) * Q + repmat( t,1,Ncorr ) ) );
            residual = dot( L, P, 1 )';
        end
        function jacobian = FJac_C_2D_LineDistance( obj, Rt )
            % Compute jacobian of error vector wrt R for observations data and certain R
            L = obj.cam_C_L;
            Q = obj.LRF_C_Q;
            
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
%         function weights  = FWeights_C_2D_LineDistance( obj, Rt )
%             residual = obj.FErr_C_2D_LineDistance( Rt );
%             residual = reshape( residual, 3, [] );
%             N = size(residual,2);
%             weights  = dot( residual, residual, 2 ) / N;
%             weights  = kron( eye(N), diag(1./weights) );
%         end
        function H = FHes_C_2D_LineDistance( obj, Rt )
            % Linear approximation to hessian in LM method
            jacobian = obj.FJac_C_2D_LineDistance( Rt );
            weights  = eye( obj.Ncorresp_C );
            H = jacobian' * weights * jacobian;
        end
        
        function [R,t] = optimizeRt_C_NonWeighted( obj )
            Fun = @(Rt) deal( obj.FErr_C_2D_LineDistance( Rt ) , obj.FJac_C_2D_LineDistance( Rt ) );
            Rt = obj.optimize( Fun, obj.Rt0, 'SE(3)', false );
            R = Rt(1:3,1:3);
            t = Rt(1:3,4);
        end
%         function [R,t] = optimizeRt_C_Weighted( obj )
%             Fun = @(Rt) deal( obj.FErr_2D_LineDistance( Rt ) ,...
%                               obj.FJac_2D_LineDistance( Rt ) ,...
%                               obj.FWeights_2D_LineDistance( Rt ) );
%             Rt = obj.optimize( Fun, obj.Rt0, 'SE(3)', true );
%             R = Rt(1:3,1:3);
%             t = Rt(1:3,4);
%         end
%         function [R,t] = optimizeRt_C_ConstWeighted( obj )
%             [R,t] = obj.optimizeRt_NonWeighted;
%             Rt0 = [R,t];
%             weights = obj.FWeights_2D_LineDistance( Rt0 );
%             
%             Fun = @(Rt) deal( obj.FErr_2D_LineDistance( Rt ) ,...
%                               obj.FJac_2D_LineDistance( Rt ) ,...
%                               weights );
%             Rt = obj.optimize( Fun, Rt0, 'SE(3)', true );
%             R = Rt(1:3,1:3);
%             t = Rt(1:3,4);
%         end
%         function [R,t] = optimizeRt_C_PreWeighted( obj )
%             [R,t] = obj.optimizeRt_NonWeighted;
%             Rt0 = [R,t];
%             
%             Fun = @(Rt) deal( obj.FErr_2D_LineDistance( Rt ) ,...
%                               obj.FJac_2D_LineDistance( Rt ) ,...
%                               obj.FWeights_2D_LineDistance( Rt ) );
%             Rt = obj.optimize( Fun, Rt0, 'SE(3)', true );
%             R = Rt(1:3,1:3);
%             t = Rt(1:3,4);
%         end
        
        function h = plotRotation_C_CostFunction( obj, R, t )
            gv  = obj.get_plot_gv( obj.plot_dist_R );
            
            FE = @(R)obj.FErr_C_2D_LineDistance( [R,t] );
            W  = eye( obj.Ncorresp_C );
            Fx = @(R,inc) expmap( inc ) * R;
            x0 = R;
            h  = obj.plotCostFunction( gv, W, FE, Fx, x0 );
        end
        
        function h = plotTranslation_C_CostFunction( obj, R, t )
            gv  = obj.get_plot_gv( obj.plot_dist_t );
            
            FE = @(t)obj.FErr_C_2D_LineDistance( [R,t] );
            W  = eye( obj.Ncorresp_C );
            Fx = @(t,inc) t + inc;
            x0 = t;
            h  = obj.plotCostFunction( gv, W, FE, Fx, x0 );
        end
        
        
        %% Get-functions
        % For translation
        function cam_L = get.cam_L( obj )
            cam_L = [obj.obs.cam_l];
            % Mask with existing measures and not-inliers
%             mask_valid = obj.mask_LRF_Q & (~obj.mask_RANSAC_t_outliers);
%             cam_L = cam_L(:, mask_valid);
            cam_L = cell2mat( cam_L );
        end
        
        function LRF_Q = get.LRF_Q( obj )
            LRF_Q = [obj.obs.LRF_q];
            % Remove non-observed data: Not necessary really (already empty)
            % Remove outliers: Set as empty observations
%             mask_valid = obj.mask_LRF_Q & (~obj.mask_RANSAC_t_outliers);
%             LRF_Q( ~mask_valid ) = [];
            % Convert from cell with empty elements to dense matrix
            LRF_Q = cell2mat( LRF_Q );
        end
        
        function cam_C_L = get.cam_C_L( obj )
            cam_C_L = reshape([obj.obs.cam_l],3,[]);
            cam_C_L = cell2mat( cam_C_L(1,:) );
        end
        
        function LRF_C_Q = get.LRF_C_Q( obj )
            LRF_C_Q = reshape([obj.obs.LRF_q],3,[]);
            LRF_C_Q = cell2mat( LRF_C_Q(1,:) );
        end

        function Rt0 = get.Rt0( obj )
            Rt0 = [ obj.R0, obj.t0 ];
        end
        function Ncorresp = get.Ncorresp( obj )
            Ncorresp = 3 * length(obj.obs);
        end
        function Ncorresp_C = get.Ncorresp_C( obj )
            Ncorresp_C = length(obj.obs);
        end
        % More functions
        function disp_N_obs( obj )
            Nobs = length(obj.obs);
            fprintf('The number of Corner complete observations is %d\n',Nobs);
        end
    end
    
end

