classdef CTrihedronOptimization < handle & CBaseOptimization
    %CTrihedronOptimization Class to store observations and later process
    %and optimize them
    %   Detailed explanation goes here
    
    properties
        % Auxiliar variables
        K   % Camera intrinsic matrix (for translation optimization)
        
        RANSAC_Rotation_threshold % Threshold for rotation error function
        RANSAC_Translation_threshold % Threshold for translation error function
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
        function obj = CTrihedronOptimization( K, RANSAC_Rotation_threshold, RANSAC_Translation_threshold,...
                debug_level,maxIters, minParamChange, minErrorChange ) % COptimization inputs
            % Set optimization parameters through parent class
            obj = obj@CBaseOptimization( debug_level, maxIters, minParamChange, minErrorChange );

            obj.obs = CTrihedronObservation.empty(1,0);
            
            if ~exist('RANSAC_Rotation_threshold','var')
                RANSAC_Rotation_threshold = 1e-1;
            end
            obj.RANSAC_Rotation_threshold = RANSAC_Rotation_threshold;
            
            if ~exist('RANSAC_Translation_threshold','var')
                RANSAC_Translation_threshold = 1e-3;
            end
            obj.RANSAC_Translation_threshold = RANSAC_Translation_threshold;
                        
            obj.K = K;
            
            % Commented because is dependent now
%             obj.mask_RANSAC_R_outliers = []; % Initialize as empty
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
        function H = FHes_Orthogonality( obj, R )
            % Linear approximation to hessian in LM method
            jacobian = obj.FJac_Orthogonality( R );
            weights  = obj.FWeights_Orthogonality( R );
            H = jacobian' * weights * jacobian;
        end
        
        function cov_data = FCov_data_Ort( obj )
            % Return covariance matrix of trihedron n's and scan v's data
            % [N|V]
            N_obs = obj.Nobs;
            CA_N = cell(1, N_obs);
            CA_V = cell(1, N_obs);
            for i=1:N_obs
                % Take observation elements which exist and are not
                % outliers
                mask_N = obj.obs(i).thereis_LRF_v & (~obj.obs(i).is_R_outlier);
                                
                % Correlated covariance of normals to planes
                A_N = obj.obs(i).cam_A_R_c_w;
                A_N = mat2cell(A_N,[3 3 3],[3 3 3]);
                A_N = A_N(mask_N,mask_N);
                CA_N{i} = cell2mat(A_N);
                
                mask_V = obj.obs(i).thereis_LRF_v & (~obj.obs(i).is_R_outlier);
                
                % Could be used non-minimal covariance
                % Input LRF_A_v is a 1x3 cell array (diagonal elements)
                A_V = obj.obs(i).LRF_A_v;
                % Remove outlier of existing part of A:V
                A_V(~mask_V) = [];
                CA_V{i} = diag(cell2mat(A_V));
            end
            A_N = blkdiag( CA_N{:} );
            A_V = blkdiag( CA_V{:} );
            A_NV = blkdiag( A_N, A_V );
            cov_data = A_NV;
        end
        function jacobian = FJac_R_NW( obj, R )
            n_cam = obj.cam_N;
            v_LRF = obj.LRF_V;
            
            if size(n_cam,2) == size(v_LRF,2) % Just to assure dimensions
                N = size(n_cam,2);
            end
            
            % Create two blocks of jacobian: J_eps_N and J_eps_V
            blockH  = cell(1,N);
            blockDN = cell(1,N);
            blockDV = cell(1,N);
            for i=1:N
                ni  = n_cam(:,i); % Trihedron normal seen from camera (3x1)
                vsi = v_LRF(:,i); % Dir vector seen from scanner (2x1)
                vci = R(:,1:2) * vsi; % Line vector seen from camera (3x1)
                
                % Phi is dC/deps = -Sum( cross(R*lsi, ni*(R*lsi)'*ni) )
                d_phi_ni  = 2 * skew(vci) * (ni*vci' + vci'*ni*eye(3));
                d_phi_vci = 2 * skew(vci) * (ni*ni') - 2 * skew(ni * vci' * ni);
                d_vci_vsi = R(:,1:2);
                d_vci_eps = - skew( R(:,1:2)*vsi );
                
                % Derivation wrt l Lie algebra instead of complete space
                d_vsi_alpha = [ -vsi(2), +vsi(1) ]';
                
                % Store blocks of H and D2
                blockH{i} = d_phi_vci * d_vci_eps;
                blockDN{i} = d_phi_ni;
                blockDV{i} = d_phi_vci * d_vci_vsi * d_vsi_alpha;
            end
            H = sum( reshape(cell2mat( blockH ),3,3,[]), 3 );
            D = [ cell2mat( blockDN ), cell2mat( blockDV ) ];
            jacobian = - H \ D;
        end
        function cov = FCov_R_NW( obj, R )
            A_NV = obj.FCov_data_Ort( );
            
            J = obj.FJac_R_NW( R );
            cov = J * A_NV * J';
        end
        function jacobian = FJac_R_W( obj, R )
            n_cam = obj.cam_N;
            v_LRF = obj.LRF_V;
            
            W = obj.FWeights_Orthogonality( R );
            
            if size(n_cam,2) == size(v_LRF,2) % Just to assure dimensions
                N = size(n_cam,2);
            end
            
            E  = obj.FErr_Orthogonality( R );
            dE_deps = obj.FJac_Orthogonality( R );
            
            % Create two blocks of jacobian: J_eps_N and J_eps_V
            blockH  = cell(1,N);
            blockDN = cell(1,N);
            blockDV = cell(1,N);
            for i=1:N
                ni  = n_cam(:,i); % Trihedron normal seen from camera (3x1)
                vsi = v_LRF(:,i); % Dir vector seen from scanner (2x1)
                vci = R(:,1:2) * vsi; % Dir vector seen from camera (3x1)
                
                % Phi is dC/deps = 2 E' * W * dE/deps
                d_ei_vci = ni';
                d_ei_ni  = vci';
                d_gi_vci = +skew( ni  );
                d_gi_ni  = -skew( vci );
                
                % Constant values related to W
                WG = ( W(i,:) * dE_deps )';
                WE = E' * W(:,i);
                d_phi_vci = 2 * ( WG * d_ei_vci + WE * d_gi_vci );
                d_phi_ni  = 2 * ( WG * d_ei_ni + WE * d_gi_ni );
                
                d_vci_vsi = R(:,1:2);
                d_vci_eps = - skew( R(:,1:2)*vsi );
                
                % Derivation wrt l Lie algebra instead of complete space
                d_vsi_alpha = [ -vsi(2), +vsi(1) ]';
                
                % Store blocks of H and D2
                blockH{i} = d_phi_vci * d_vci_eps;
                blockDN{i} = d_phi_ni;
                blockDV{i} = d_phi_vci * d_vci_vsi * d_vsi_alpha;
            end
            H = sum( reshape(cell2mat( blockH ),3,3,[]), 3 );
            D = [ cell2mat( blockDN ), cell2mat( blockDV ) ];
            jacobian = - H \ D;
        end
        function cov = FCov_R_W( obj, R )
            A_NV = obj.FCov_data_Ort( );
            
            J = obj.FJac_R_W( R );
            cov = J * A_NV * J';
        end
        
        R = optimizeRotation_NonWeighted( obj )
        R = optimizeRotation_Weighted( obj )
        R = optimizeRotation_DiagWeighted( obj )
        J = FJac_optimizeRotation_NonWeighted( obj, R )
        
        % TODO: Implement optimizeRotation_Covariance from optimRotation.m
        function h = plotRotationCostFunction( obj, R )
            gv  = obj.get_plot_gv( obj.plot_dist_R );
            
            FE = @(R)obj.FErr_Orthogonality( R );
            W  = obj.FWeights_Orthogonality( R );
            Fx = @(R,inc) expmap( inc ) * R;
            x0 = R;
            h  = obj.plotCostFunction( gv, W, FE, Fx, x0 );
        end
%         h = plotRotationCostFunction( obj, R )
        
        plot( obj )
        
        % For Translation
        % 3D error functions
        function residual = FErr_3D_PlaneDistance( obj, R, t )
            % Compute error vector for observations data, R and t
            % N - (3x...) 3D normals to reprojection planes from camera center
            % through image lines
            % Q - (2x...) 2D LRF intersection points
            N = obj.cam_reprN;
            Q = obj.LRF_Q;
            
            s = size(N,2); % Number of correspondences
            
            residual = dot( N, R(:,1:2)*Q + repmat(t,1,s), 1)';
        end
        function jacobian = FJac_3D_PlaneDistance( obj, R, t )
            % Compute jacobian of error vector wrt R for observations data and certain R
            N = obj.cam_reprN;
            jacobian = N';
        end
        weights = FWeights_3D_PlaneDistance( obj, R, t )
        function H = FHes_3D_PlaneDistance( obj, R, t )
            % Linear approximation to hessian in LM method
            jacobian = obj.FJac_3D_PlaneDistance( R, t );
            weights  = obj.FWeights_3D_PlaneDistance( R, t );
            H = jacobian' * weights * jacobian;
        end
        function cov_data = FCov_data_3D( obj )
            N_obs = obj.Nobs;
            CA_N = cell(1, N_obs);
            CA_Q = cell(1, N_obs);
            
            for i=1:N_obs
                % Take observation elements which exist and are not
                % outliers
                mask_N = obj.obs(i).thereis_LRF_q & (~obj.obs(i).is_t_outlier);
                               
                % Correlated covariance of normals to planes
                A_N = obj.obs(i).cam_A_reprN;
                A_N = mat2cell(A_N,[3 3 3],[3 3 3]);
                A_N = A_N(mask_N,mask_N);
                CA_N{i} = cell2mat(A_N);
                
                mask_Q = obj.obs(i).thereis_LRF_q & (~obj.obs(i).is_t_outlier);
                % Could be used non-minimal covariance
                % Input LRF_A_v is a 1x3 cell array (diagonal elements)
                A_Q = obj.obs(i).LRF_A_q;
                % Remove outlier of existing part of A:V
                A_Q(~mask_Q) = [];
                CA_Q{i} = blkdiag( A_Q{:} );
            end
            A_N = blkdiag( CA_N{:} );
            A_Q = blkdiag( CA_Q{:} );
            A_NQ = blkdiag( A_N, A_Q );
            cov_data = A_NQ;
        end
        function jacobian = FJac_t_3D_W( obj, R, t )
            N = obj.cam_reprN;
            Q = obj.LRF_Q;
            
            W = obj.FWeights_3D_PlaneDistance( R, t );
            
            if size(N,2) == size(Q,2) % Just to assure dimensions
                Ncorresp = size(N,2);
            end
            
            E  = obj.FErr_3D_PlaneDistance( R, t );
            dE_dt = obj.FJac_3D_PlaneDistance( R, t );
            
            % Create two blocks of jacobian: J_eps_N and J_eps_Q
            blockH  = cell(1,Ncorresp);
            blockDN = cell(1,Ncorresp);
            blockDQ = cell(1,Ncorresp);
            for i=1:Ncorresp
                ni  = N(:,i); % Back-projected plane normal seen from camera (3x1)
                qsi = Q(:,i); % Corner point seen from scanner (2x1)
                qci = R(:,1:2) * qsi + t; % Corner point seen from camera (3x1)
                
                % Phi is dC/deps = 2 E' * W * dE/deps
                d_ei_qci = ni';
                d_ei_ni  = qci';
                d_gi_qci = zeros(3,3); % Non-dependent
                d_gi_ni  = eye(3);
                
                % Constant values related to W
                WG = ( W(i,:) * dE_dt )';
                WE = E' * W(:,i);
                d_phi_qci = 2 * ( WG * d_ei_qci + WE * d_gi_qci );
                d_phi_ni  = 2 * ( WG * d_ei_ni + WE * d_gi_ni );
                
                d_qci_qsi = R(:,1:2);
                d_qci_eps = - skew( R(:,1:2)*qsi );
                
                % Store blocks of H and D2
                blockH{i} = d_phi_qci * d_qci_eps;
                blockDN{i} = d_phi_ni;
                blockDQ{i} = d_phi_qci * d_qci_qsi;
            end
            H = sum( reshape(cell2mat( blockH ),3,3,[]), 3 );
            D = [ cell2mat( blockDN ), cell2mat( blockDQ ) ];
            jacobian = - H \ D;
        end
        function cov = FCov_t_3D_W( obj, R, t )
            A_NQ = obj.FCov_data_3D;
            
            J = obj.FJac_t_3D_W( R, t );
            cov = J * A_NQ * J';
        end
        
        function h = plotTranslation_3D_CostFunction( obj, R, t )
            gv  = obj.get_plot_gv( obj.plot_dist_t );
            
            FE = @(t)obj.FErr_3D_PlaneDistance( R, t );
            W  = obj.FWeights_3D_PlaneDistance( R, t );
            Fx = @(t,inc) t + inc;
            x0 = t;
            h  = obj.plotCostFunction( gv, W, FE, Fx, x0 );
        end
        
        % 2D error functions
        function residual = FErr_2D_LineDistance( obj, R, t )
            % Compute error vector for observations data and certain Rt
            L = obj.cam_L;
            Q = obj.LRF_Q;
            
            % Convert L to pixel space and normalize line:
            L = obj.K' \ L;
            L = L ./ repmat( sqrt(sum(L(1:2,:).^2,1)), 3,1 );

            Ncorr = size( Q,2 );
            P = hnormalise( obj.K * ( R(:,1:2) * Q + repmat( t,1,Ncorr ) ) );
            residual = dot( L, P, 1 )';
        end
        function jacobian = FJac_2D_LineDistance( obj, R, t )
            % Compute jacobian of error vector wrt R for observations data and certain R
            L = obj.cam_L;
            Q = obj.LRF_Q;
            
            % Convert L to pixel space and normalize line:
            L = obj.K' \ L;
            L = L ./ repmat( sqrt(sum(L(1:2,:).^2,1)), 3,1 );
                              
            Ncorr = size( Q,2 );
            jacobian = zeros(Ncorr,3);
            for i=1:Ncorr
                % Homogeneous 3D points in camera frame (pixel units)
                l = L(:,i);
                q = Q(:,i);
                ph = obj.K * ( R(:,1:2) * q + t );
                J_hnormalise = 1/ph(3)^2 * [ ph(3) 0 -ph(1)
                                             0 ph(3) -ph(2)
                                               zeros(1,3)   ];
                J_ph_t = obj.K;
                
                jacobian(i,:) = l' * J_hnormalise * J_ph_t;
            end
        end
        weights = FWeights_2D_LineDistance( obj, R, t )
        function H = FHes_2D_LineDistance( obj, R, t )
            % Linear approximation to hessian in LM method
            jacobian = obj.FJac_2D_LineDistance( R, t );
            weights  = obj.FWeights_2D_LineDistance( R, t );
            H = jacobian' * weights * jacobian;
        end
        
        function h = plotTranslation_2D_CostFunction( obj, R, t )
            gv  = obj.get_plot_gv( obj.plot_dist_t );
            
            FE = @(t)obj.FErr_2D_LineDistance( R, t );
            W  = obj.FWeights_2D_LineDistance( R, t );
            Fx = @(t,inc) t + inc;
            x0 = t;
            h  = obj.plotCostFunction( gv, W, FE, Fx, x0 );
        end
        
        function t = optimizeTranslation_2D_NonWeighted( obj, R )
            Fun = @(t) deal( obj.FErr_2D_LineDistance(R,t) , obj.FJac_2D_LineDistance(R,t) );
            t = obj.optimize( Fun, obj.t0, 'Rn', false );
        end
        function t = optimizeTranslation_2D_Weighted( obj, R )
            Fun = @(t) deal( obj.FErr_2D_LineDistance(R,t) ,...
                             obj.FJac_2D_LineDistance(R,t) ,...
                             obj.FWeights_2D_LineDistance(R,t) );
            t = obj.optimize( Fun, obj.t0, 'Rn', true );
        end

        % For global optimization (R+t)
        % Orthogonality and 3D error functions
        function residual = FErr_Global_Ort_3D( obj, R, t )
            % Compute error vector for observations data, R and t
            % N - (3x...) 3D normals to reprojection planes from camera center
            % through image lines
            % Q - (2x...) 2D LRF intersection points
            res_ort = obj.FErr_Orthogonality( R );
            res_3D  = obj.FErr_3D_PlaneDistance( R, t );
            residual = [ res_ort ; res_3D ];
        end
        function jacobian = FJac_Global_Ort_3D( obj, R, t )
            % Compute jacobian of error vector wrt R for observations data and certain R
            % Jac e_ort wrt R and t
            jac_ort_R = obj.FJac_Orthogonality( R );
            jac_ort_t = zeros( size(jac_ort_R,1), 3 );
            % Jac e_3D wrt R and t
            N = obj.cam_reprN;
            Q = obj.LRF_Q;
            jac_3D_R = cross( R(1:3,1:2) * Q, N, 1 )';
            jac_3D_t = N';
            jacobian = [ jac_ort_R , jac_ort_t ;
                         jac_3D_R  , jac_3D_t  ];
        end
        function weights = FWeights_Global_Ort_3D( obj, R, t )
            W_ort = obj.FWeights_Orthogonality( R );
            W_3D  = obj.FWeights_3D_PlaneDistance( R, t );
            weights = blkdiag( W_ort, W_3D );
        end
        function H = FHes_Global_Ort_3D( obj, R, t )
            % Linear approximation to hessian in LM method
            jacobian = obj.FJac_Global_Ort_3D( R, t );
            weights  = obj.FWeights_Global_Ort_3D( R, t );
            H = jacobian' * weights * jacobian;
        end
        function jacobian = FJac_Rt_Global_Ort_3D( obj, R, t )
            jac_R_ort = obj.FJac_R_W( R );
            jac_t_ort = zeros( 3, size(jac_R_ort,2) );
            jac_t_3D  = obj.FJac_t_3D_W( R, t );
            
            N = obj.cam_reprN;
            Q = obj.LRF_Q;
            
            W = obj.FWeights_3D_PlaneDistance( R, t );
            
            if size(N,2) == size(Q,2) % Just to assure dimensions
                Ncorresp = size(N,2);
            end
            
            E  = obj.FErr_3D_PlaneDistance( R, t );
            
            % A little risky: not completely checked
            dE_deps = cross( R(1:3,1:2) * Q, N, 1 )';
            
            % Create two blocks of jacobian: J_eps_N and J_eps_Q
            blockH  = cell(1,Ncorresp);
            blockDN = cell(1,Ncorresp);
            blockDQ = cell(1,Ncorresp);
            for i=1:Ncorresp
                ni  = N(:,i); % Back-projected plane normal seen from camera (3x1)
                qsi = Q(:,i); % Corner point seen from scanner (2x1)
                qci = R(:,1:2) * qsi + t; % Corner point seen from camera (3x1)
                
                % Phi is dC/deps = 2 E' * W * dE/deps
                d_ei_qci = ni';
                d_ei_ni  = qci';
                d_gi_qci = +skew(ni); % Non-dependent
                d_gi_ni  = +skew(R(:,1:2)*qsi)';
                
                % Constant values related to W
                WG = ( W(i,:) * dE_deps )';
                WE = E' * W(:,i);
                d_phi_qci = 2 * ( WG * d_ei_qci + WE * d_gi_qci );
                d_phi_ni  = 2 * ( WG * d_ei_ni + WE * d_gi_ni );
                
                d_qci_qsi = R(:,1:2);
                d_qci_eps = - skew( R(:,1:2)*qsi );
                
                % Store blocks of H and D2
                blockH{i} = d_phi_qci * d_qci_eps;
                blockDN{i} = d_phi_ni;
                blockDQ{i} = d_phi_qci * d_qci_qsi;
            end
            H = sum( reshape(cell2mat( blockH ),3,3,[]), 3 );
            D = [ cell2mat( blockDN ), cell2mat( blockDQ ) ];
            jac_R_3D = - H \ D;
            
            jacobian = [ jac_R_ort , jac_R_3D ;
                         jac_t_ort , jac_t_3D ];
        end
        function cov = FCov_Global_Ort_3D( obj, R, t )
            A_NV = obj.FCov_data_Ort;
            A_NQ = obj.FCov_data_3D;
            A_NVNQ = blkdiag( A_NV, A_NQ );
            
            J = obj.FJac_Rt_Global_Ort_3D( R, t );
            cov = J * A_NVNQ * J';
        end
        
        function h = plotTranslation_Global_CostFunction( obj, R, t )
            gv  = obj.get_plot_gv( obj.plot_dist_t );
            
            FE = @(t)obj.FErr_Global_Ort_3D( R, t );
            W  = obj.FWeights_Global_Ort_3D( R, t );
            Fx = @(t,inc) t + inc;
            x0 = t;
            h  = obj.plotCostFunction( gv, W, FE, Fx, x0 );
        end
        function h = plotRotation_Global_CostFunction( obj, R, t, dist )
            if ~exist( 'dist','var' )
                dist = obj.plot_dist_R;
            end
            gv  = obj.get_plot_gv( dist );
            
            FE = @(R)obj.FErr_Global_Ort_3D( R, t );
            W  = obj.FWeights_Global_Ort_3D( R, t );
            Fx = @(R,inc) expmap( inc ) * R;
            x0 = R;
            h  = obj.plotCostFunction( gv, W, FE, Fx, x0 );
        end
        
        function [R,t] = optimizeGlobal_Ort_3D( obj, R, t )
            Fun = @(Rt) deal( obj.FErr_Global_Ort_3D(Rt(1:3,1:3),Rt(1:3,4)) ,...
                             obj.FJac_Global_Ort_3D(Rt(1:3,1:3),Rt(1:3,4)) ,...
                             obj.FWeights_Global_Ort_3D(Rt(1:3,1:3),Rt(1:3,4)) );
%             Rt = obj.optimize( Fun, [obj.R0, obj.t0], 'SE(3)', true );
            Rt = obj.optimize( Fun, [R, t], 'SE(3)', true );
            R = Rt(1:3,1:3); t = Rt(1:3,4);
        end
        
        % For global optimization (R+t)
        % 3D back-projected plane distance as error function (without
        % Orthogonality)
        function residual = FErr_Global_3D( obj, R, t )
            % Compute error vector for observations data, R and t
            % N - (3x...) 3D normals to reprojection planes from camera center
            % through image lines
            % Q - (2x...) 2D LRF intersection points
            res_3D  = obj.FErr_3D_PlaneDistance( R, t );
            residual = res_3D;
        end
        function jacobian = FJac_Global_3D( obj, R, t )
            % Compute jacobian of error vector wrt R for observations data and certain R
            % Jac e_3D wrt R and t
            N = obj.cam_reprN;
            Q = obj.LRF_Q;
            jac_3D_R = cross( R(1:3,1:2) * Q, N, 1 )';
            jac_3D_t = N';
            jacobian = [ jac_3D_R  , jac_3D_t  ];
        end
        function weights = FWeights_Global_3D( obj, R, t )
            W_3D  = obj.FWeights_3D_PlaneDistance( R, t );
            weights = W_3D;
        end
        function H = FHes_Global_3D( obj, R, t )
            % Linear approximation to hessian in LM method
            jacobian = obj.FJac_Global_3D( R, t );
            weights  = obj.FWeights_Global_3D( R, t );
            H = jacobian' * weights * jacobian;
        end
        function jacobian = FJac_Rt_Global_3D( obj, R, t )
            % TOCHECK!!! Copied from Ort+3D
            jac_t_3D  = obj.FJac_t_3D_W( R, t );
            
            N = obj.cam_reprN;
            Q = obj.LRF_Q;
            
            W = obj.FWeights_3D_PlaneDistance( R, t );
            
            if size(N,2) == size(Q,2) % Just to assure dimensions
                Ncorresp = size(N,2);
            end
            
            E  = obj.FErr_3D_PlaneDistance( R, t );
            
            % A little risky: not completely checked
            dE_deps = cross( R(1:3,1:2) * Q, N, 1 )';
            
            % Create two blocks of jacobian: J_eps_N and J_eps_Q
            blockH  = cell(1,Ncorresp);
            blockDN = cell(1,Ncorresp);
            blockDQ = cell(1,Ncorresp);
            for i=1:Ncorresp
                ni  = N(:,i); % Back-projected plane normal seen from camera (3x1)
                qsi = Q(:,i); % Corner point seen from scanner (2x1)
                qci = R(:,1:2) * qsi + t; % Corner point seen from camera (3x1)
                
                % Phi is dC/deps = 2 E' * W * dE/deps
                d_ei_qci = ni';
                d_ei_ni  = qci';
                d_gi_qci = +skew(ni); % Non-dependent
                d_gi_ni  = +skew(R(:,1:2)*qsi)';
                
                % Constant values related to W
                WG = ( W(i,:) * dE_deps )';
                WE = E' * W(:,i);
                d_phi_qci = 2 * ( WG * d_ei_qci + WE * d_gi_qci );
                d_phi_ni  = 2 * ( WG * d_ei_ni + WE * d_gi_ni );
                
                d_qci_qsi = R(:,1:2);
                d_qci_eps = - skew( R(:,1:2)*qsi );
                
                % Store blocks of H and D2
                blockH{i} = d_phi_qci * d_qci_eps;
                blockDN{i} = d_phi_ni;
                blockDQ{i} = d_phi_qci * d_qci_qsi;
            end
            H = sum( reshape(cell2mat( blockH ),3,3,[]), 3 );
            D = [ cell2mat( blockDN ), cell2mat( blockDQ ) ];
            jac_R_3D = - H \ D;
            
            jacobian = [ jac_R_ort , jac_R_3D ;
                         jac_t_ort , jac_t_3D ];
        end
        function cov = FCov_Global_3D( obj, R, t )
            A_NQ = obj.FCov_data_3D;
            
            J = obj.FJac_Rt_Global_3D( R, t );
            cov = J * A_NV * J';
        end
        
        function [R,t] = optimizeGlobal_3D( obj, R, t )
            Fun = @(Rt) deal( obj.FErr_Global_3D(Rt(1:3,1:3),Rt(1:3,4)) ,...
                             obj.FJac_Global_3D(Rt(1:3,1:3),Rt(1:3,4)) ,...
                             obj.FWeights_Global_3D(Rt(1:3,1:3),Rt(1:3,4)) );
%             Rt = obj.optimize( Fun, [obj.R0, obj.t0], 'SE(3)', true );
            Rt = obj.optimize( Fun, [R, t], 'SE(3)', true );
            R = Rt(1:3,1:3); t = Rt(1:3,4);
        end
        
        
%         t = optimizeTranslation_2D_NonWeighted( obj, R )
%         t = optimizeTranslation_2D_Weighted( obj, R )
        t = optimizeTranslation_3D_NonWeighted( obj, R )
        t = optimizeTranslation_3D_Weighted( obj, R )
        
        function showSamplingSphere( obj )
            N = obj.cam_N;
            L = obj.LRF_V;
            
            vecL = zeros(3, size(N,2));
            for i=1:size(N,2) % For representing measures in 3D space
                n = N(:,i);
                %     J_n_y = Householder(n) * [eye(2) , [0 0]']';
                J_n_y = Householder(n);
                J_n_y = J_n_y(:,1:2);
                vecL(:,i) = J_n_y * L(:,i);
            end
            figure('Name','Sampling map'), hold on
            scatter3( N(1,:), N(2,:), N(3,:), 'r' )
            quiver3( N(1,:), N(2,:), N(3,:), vecL(1,:), vecL(2,:), vecL(3,:), 0.1, 'b' )
            axis equal, rotate3d on
        end
        
        %% Get-functions
        % For rotation
        function cam_N = get.cam_N( obj )
            cam_N = [obj.obs(1:obj.Nobs).cam_R_c_w];
            % Mask with existing measures and not-inliers
            mask_valid = obj.mask_LRF_V & (~obj.mask_RANSAC_R_outliers);
            cam_N = cam_N(:, mask_valid);
        end        
        function LRF_V = get.LRF_V( obj )
            LRF_V = [obj.obs(1:obj.Nobs).LRF_v];
            % Remove non-observed data: Not necessary really (already empty)
            % Remove outliers: Set as empty observations
            mask_valid = obj.mask_LRF_V & (~obj.mask_RANSAC_R_outliers);
            LRF_V( ~mask_valid ) = [];
            % Convert from cell with empty elements to dense matrix
            LRF_V = cell2mat( LRF_V );
        end
        % For translation
        function cam_L = get.cam_L( obj )
            cam_L = [obj.obs(1:obj.Nobs).cam_l];
            % Mask with existing measures and not-inliers
            mask_valid = obj.mask_LRF_Q & (~obj.mask_RANSAC_t_outliers);
            cam_L = cam_L(:, mask_valid);
        end
        function cam_reprN = get.cam_reprN( obj )
            cam_reprN = [obj.obs(1:obj.Nobs).cam_reprN];
            % Mask with existing measures and not-inliers
            mask_valid = obj.mask_LRF_Q & (~obj.mask_RANSAC_t_outliers);
            cam_reprN = cam_reprN(:, mask_valid);
        end
        function LRF_Q = get.LRF_Q( obj )
            LRF_Q = [obj.obs(1:obj.Nobs).LRF_q];
            % Remove non-observed data: Not necessary really (already empty)
            % Remove outliers: Set as empty observations
            mask_valid = obj.mask_LRF_Q & (~obj.mask_RANSAC_t_outliers);
            LRF_Q( ~mask_valid ) = [];
            % Convert from cell with empty elements to dense matrix
            LRF_Q = cell2mat( LRF_Q );
        end
        function mask_LRF_V = get.mask_LRF_V( obj )
            mask_LRF_V = [obj.obs(1:obj.Nobs).thereis_LRF_v];
        end
        function mask_LRF_Q = get.mask_LRF_Q( obj )
            mask_LRF_Q = [obj.obs(1:obj.Nobs).thereis_LRF_q];
        end
        function mask_RANSAC_R_outliers = get.mask_RANSAC_R_outliers( obj )
            mask_RANSAC_R_outliers = [obj.obs(1:obj.Nobs).is_R_outlier];
        end
        function mask_RANSAC_t_outliers = get.mask_RANSAC_t_outliers( obj )
            mask_RANSAC_t_outliers = [obj.obs(1:obj.Nobs).is_t_outlier];
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

