classdef CCheckerboardOptimization < handle & CBaseOptimization
    %CCornerOptimization Class to store observations and later process
    %and optimize them
    %   Detailed explanation goes here
    
    properties       
        % Auxiliar variables
        K   % Camera intrinsic matrix (for translation optimization)
    end
    
    properties (SetAccess=private)
        % Optimized variables
        R
        t
        
        inliers     % Field to store inliers found by Vasconcelos to use them in Zhang
    end
    
    properties (Dependent) % Array values collected from set of observations
        plane_T     % Grid poses in 3D (TODO: wrt Cam? Cam to Plane?)
        LRF_pts     % LRF points intersecting plane (in [mm]!!!)
    end
    
    methods
        %% Constructor
        function obj = CCheckerboardOptimization( K, ...
                debug_level,maxIters, minParamChange, minErrorChange ) % COptimization inputs
            % Set optimization parameters through parent class
            obj = obj@CBaseOptimization( debug_level, maxIters, minParamChange, minErrorChange );

            obj.obs = CCheckerboardObservation.empty(1,0);
            
            obj.K = K;
        end
                        
        %% Compute initial estimate
        function obj = setInitialRotation( obj, R )
            obj.R = R;
        end
        
        function obj = setInitialTranslation( obj, t )
            obj.t = t;
        end
        
        %% Optimization functions and methods
        function [R,t] = optimizeRt_Vasc(obj)
            Ti = [obj.R, obj.t];
            % Convert from Cam frame to LRF frame (t in [mm])
            Ti = [Ti(1:3,1:3)', -1000*Ti(1:3,1:3)'*Ti(1:3,4)];
            [T, ~, obj.inliers,~,~] = lccMinSol(obj.plane_T,obj.LRF_pts,Ti); % Added initial estimate to save time
            % * Store inliers (planes) to speed Zhang
            
            % The criterion is different for Vasconcelos (T_S_C) -> Convert
            T = inv(T);
            R = T(1:3,1:3);
            t = T(1:3,4) / 1000; % Change mm to m
        end
        
        function [R,t] = optimizeRt_Zhang(obj)
            % Vasconcelos has to be applied first to obtain inliers
            if isempty(obj.inliers)
                error('Zhang: No inliers. Vasconcelos has to be applied fist to get inliers');
            end
            Ti = [obj.R, obj.t];
            % Convert from Cam frame to LRF frame (t in [mm])
            Ti = [Ti(1:3,1:3)', -1000*Ti(1:3,1:3)'*Ti(1:3,4)];
            [T, ~,~,~,~] = lccZhang(obj.plane_T,obj.LRF_pts, Ti,obj.inliers);
            obj.inliers = NaN; % Set NaN to check Zhang used them
            % The criterion is different for Zhang (T_S_C) -> Convert
            T = inv(T);
            R = T(1:3,1:3);
            t = T(1:3,4) / 1000; % Change mm to m
        end
        
        
        %% Get-functions
        function plane_T = get.plane_T( obj )
            plane_T = [obj.obs(1:obj.Nobs).plane_T];
            plane_T = reshape( plane_T, 4, 4, [] );
        end
        function LRF_pts = get.LRF_pts( obj )
            LRF_pts = {obj.obs(1:obj.Nobs).LRF_pts};
        end
        
        % More functions
        function disp_N_obs( obj )
%             Nobs = length(obj.obs);
            Nobs = obj.Nobs;
            fprintf('The number of Checkerboard complete observations is %d\n',Nobs);
        end
    end
    
end