classdef CTrihedronSolver < handle
    %CTrihedronSolver Class for solving trihedron problem
    %   Detailed explanation goes here
    
    properties
        % Input
        Nbp % Back-projected planes normals
        A_Nbp % Back-projected planes covariance matrices
        K   % Camera intrinsic calibration
        c   % Cos of angles among trihedron directions
        
        % Intermediate variables
        rho
        lam
        Om
        U
        OmU
        B
        D
        
        % Output
        V
        Am_V
        obj_V
        
        % Solution selection variables
        Dproj
        v_im
    end
    
    methods
        function this = CTrihedronSolver( Nbp, K, c )
            this.Nbp = Nbp;
            this.K   = K;
            if ~exist('c','var')
                c = zeros(3,1);
            end
            this.c = c;
            
            % Intermediate variables
            this.B = cell(1,3);
            
            % Auxiliar variables
            this.Dproj = cell(1,3);
            this.v_im  = cell(1,3);
        end
        
        % Auxiliar loading functions
        function loadXi( this, xi )
            this.Dproj = deal( DprojectionP2( xi.c.X ) );
            this.v_im  = { xi.v(1).X , xi.v(2).X , xi.v(3).X };
        end
        
        function loadSegments( this, segments )
            for k=1:3
                this.Dproj{k} = DprojectionP2( segments(k).p1 );
                this.v_im{k}  = segments(k).v;
            end
        end
        
        function loadCovariance( this, C_A_N )
            this.A_Nbp = C_A_N;
        end
        
        % Auxiliar reading functions
        function obj_V = getManifoldObject( this )
            obj_V = this.obj_V;
        end
        
        % Uncertainty estimation
        function computeCovariance( this )
            % Compute covariance by 1st propagation from
            % normals of interpretation planes
            if ~any(this.c)
                this.obj_V = Manifold.SO3( this.V );
                
                    % Compute covariance of trihedron normals
                    J_Phi_eps = cross( this.V, this.Nbp, 1 )';
                    J_Phi_Nbp = mat2cell( this.V, 3, [1 1 1] );
                    J_Phi_Nbp = blkdiag( J_Phi_Nbp{:} )';
                    J_eps_Nbp = - J_Phi_eps \ J_Phi_Nbp;
                    
                    if numel( this.A_Nbp ) == 3
                        % Diagonal non-correlated case
                        A_N = blkdiag( this.A_Nbp{:} );
                    elseif all( size( this.A_Nbp ) == [3 3] )
                        % Full (correlated) case
                        A_N = cell2mat( this.A_Nbp );
                    else
                        error('Wrong A_Nbp property, check');
                    end
                    this.obj_V.setMinimalCov( ...
                        J_eps_Nbp * A_N * J_eps_Nbp' );
                    this.Am_V = this.obj_V.A_x;
            end
        end
        
        % Solver functions
        function checkJunction( this )           
            % Check if 3-junction
            if abs( det( this.Nbp ) ) < 1e-10
                % Consider simple meeting lines
                this.getReducedJunction;
                this.chooseTrihedronNormals;
                return
            end
            
            % Compute intersection direction t (of planes 1 and 2)
            n = mat2cell(this.Nbp,3,[1 1 1]);
            t = snormalize( cross( n{1}, n{2} ) );
            Om_t = null(t');
            
            % delta is the reduced vector for n3 direction,
            % n3 = null(t') * delta
            % Define new bilinear form
            nm = n{3};
            M = (nm' * t)^2 / (n{1}' * skew(t)*skew(t) * n{2})...
                * ( skew(t) * n{2} * n{1}' * skew(t) ) + nm*nm';
            M = Om_t' * skew(t) * M * skew(t) * Om_t;
            M = 0.5 * (M+M'); % Symmetrize
            
            [U,D] = eig( M ); % Diagonalize matrix
            if sign(D(1,1)) == sign(D(2,2))
                warning('Not different signs');
                keyboard
            end
            n3 = Om_t * U * snormalize(...
                  [ +sqrt(abs(D(2,2))), +sqrt(abs(D(1,1)))
                    -sqrt(abs(D(2,2))), +sqrt(abs(D(1,1))) ]');
            % Check good solution (for each a0 sign)
            a12 = n{1}'*skew(t)*skew(t)*n{2};
            a23 = n{2}'*skew(t)*skew(t)*n3;
            a31 = n{1}'*skew(t)*skew(t)*n3;
            a0  = sqrt( [a12 a12] .* a23 .* a31 );
            err1 = +nm'*t*a0 + a12 * nm'*skew(t)*n3;
            err2 = -nm'*t*a0 + a12 * nm'*skew(t)*n3;
            [~,ind(1)] = min(abs(err1));
            [~,ind(2)] = min(abs(err2));
            
            a0_sign = [+1 -1];
            for ii = 1:2
                this.Nbp(:,3) = n3(:,ind(ii));
                
                this.getReducedJunction;
                Det = this.setSigns( diag([a0_sign(ii) 1]) * this.rho );
                if Det == 1
                    break
                end
            end

        end

        function [rho, Om] = getReducedJunction( this )
            % [rho, Om, D] = CTrihedronSolver.getReducedJunction
            % Solve the system of equations for trihedron directions
            % in reduced bases
            
            if abs(det( this.Nbp )) > 1e-10
                error('Use getReducedJunction only with meeting corners');
            end
            
            % Obtain null spaces (Om) in which v lies (orthogonality condition)
            n = mat2cell(this.Nbp,3,[1 1 1]); % Reshape to cell array
            t = snormalize( cross( n{1}, n{2} ) );
            Om = cell(1,3);
            Om{1} = [ t, cross(t,n{1}) ];
            Om{2} = [ t, cross(t,n{2}) ];
            Om{3} = [ t, cross(t,n{3}) ];
            
            % Compute bilinear forms matrices
            % Bilinear forms become diagonal for chosen bases, with form
            % D_ij = [ 1 0 ; 0 -alpha_ij ];
            d = zeros(1,3);
            for k=1:3
                ij = setdiff(1:3,k); i=ij(1); j=ij(2);
                d(k) = n{i}'*skew(t)*skew(t)*n{j};
            end            
            d_m = sqrt( prod(d) );
            rho = snormalize( [ d_m d_m d_m ; d ] );
            
            % Store class intermediate results
            this.rho = rho;
            this.Om  = Om;
        end
        
        function [rho, OmU, D] = getReducedParams( this )
            % [rho, OmU, D] = CTrihedronSolver.solveReducedParams
            % Solve the system of equations for trihedron directions
            % in reduced bases
            
            % Obtain null spaces (Om) in which v lies (orthogonality condition)
            cell_Nbp = mat2cell(this.Nbp,3,[1 1 1]); % Reshape to cell array
            Om = cell(1,3);
            for k=1:3
                Om{k} = null( cell_Nbp{k}' );
            end
            
            % Compute bilinear forms matrices
            B = cell(3,3);
            for k=1:3
                ij = setdiff(1:3,k); i=ij(1); j=ij(2);
                B{i,j} = Om{i}' * Om{j};
                B{j,i} = Om{j}' * Om{i};
            end
            
            % Store bilinear forms in object
            this.B{3} = B{1,2};
            this.B{1} = B{2,3};
            this.B{2} = B{3,1};
            
            % Diagonalize bilinear forms
            U = cell(1,3);
            [U{1},S{3},V{2}] = svd( B{1,2} );
            [U{2},S{1},V{3}] = svd( B{2,3} );
            [U{3},S{2},V{1}] = svd( B{3,1} );
                        
            % Make each base dextrorotatory
            for k=1:3
                U{k} = U{k} * diag( [1 det(U{k})] );
            end
            
            D = cell(1,3);
            d = zeros(1,3);
            for k=1:3
                ij = setdiff(1:3,k); i=ij(1); j=ij(2);
                D{k} = U{i}'*B{i,j}*U{j};
                
                signD11 = sign(D{k}(1,1));
                d(k) = -D{k}(2,2) * signD11;
            end
            
            d_m = sqrt( prod(d) );
            rho = snormalize( [ d_m d_m d_m ; d ] );
            
            % Return basis changes
            OmU = cell(1,3);
            for k=1:3
                OmU{k} = Om{k} * U{k};
            end
            
            % Store class intermediate results
            this.rho = rho;
            this.Om  = Om;
            this.U   = U;
            this.OmU = OmU;
            this.D   = D;
        end
        
        function rho = chooseTrihedronNormals( this )
            % CTrihedronSolver.chooseTrihedronNormales
            % Using information about segments sign in image,
            % choose dextrorotatory solution
            
            % Two possible signs solutions:
            rho_{1} = this.rho;
            rho_{2} = diag([-1 1]) * this.rho;
            
            % Check projection sign and rotation determinant
            for k=1:2
                rho = rho_{k};
                n = cell(1,3);
                for i=1:3
                    s = sign( this.v_im{i}' *...
                              this.Dproj{i} *...
                              this.K * this.Om{i} *...
                              rho(:,i) );
                    rho(:,i) = s * rho(:,i);
                    n{i} = this.Om{i} * rho(:,i);
                end
                R_tri = cell2mat( n );
                
                if abs( det(R_tri) - 1 ) < 0.5
                    break % Exit loop with current rho output value
                end
            end
            
            this.rho = rho; % Store good solution
        end
        
        function determinant = setSigns( this, rho )
            % CTrihedronSolver.setSigns
            % Using information about segments sign in image,
            % choose dextrorotatory solution
            
            n = cell(1,3);
            for i=1:3
                s = sign( this.v_im{i}' *...
                    this.Dproj{i} *...
                    this.K * this.Om{i} *...
                    rho(:,i) );
                rho(:,i) = s * rho(:,i);
                n{i} = this.Om{i} * rho(:,i);
            end
            R_tri = cell2mat( n );
            determinant = sign( det(R_tri) );
            
            this.rho = rho; % Store new solution
        end
        
        function rho = refineReducedSystem( ~, B, c, rho0 )
            
            % Definition of necessary functions
            null2 = @(x)[-x(2); x(1)]; % Compute orthogonal vector in positive orientation
            function [F,JF] = Fun( rho )
                rho = mat2cell( rho, 2,[1 1 1] );
                F  = [ rho{1}'*B{3}*rho{2} - c(3)
                       rho{2}'*B{1}*rho{3} - c(1)
                       rho{3}'*B{2}*rho{1} - c(2) ];
                JF = [ rho{2}'*B{3}*null2(rho{1}') , rho{1}'*B{3}*null2(rho{2}') , 0
                       0 , rho{3}'*B{1}*null2(rho{2}') , rho{2}'*B{1}*null2(rho{3}')
                       rho{3}'*B{2}*null2(rho{1}') , 0 , rho{1}'*B{2}*null2(rho{3}') ];
            end
            
            obj_rho = Manifold.Dyn( Manifold.S1(rho0{1}),...
                                    Manifold.S1(rho0{2}),...
                                    Manifold.S1(rho0{3}) );
            
            [ obj_rho, Err ] = Newton_solve( @(var)Fun(reshape(var.X,2,3)),...
                obj_rho, 'space', 'Manifold', 'debug', false );
            rho = reshape( obj_rho.X, 2,3 );
            rho = mat2cell( rho, 2,[1 1 1] );
        end
        
        function rho2lam( this )
            % rho2lam( this )
            % Undo simultaneous diagonalization
            for k=1:3
                this.lam{k} = this.U{k} * this.rho(:,k);
            end
        end
        
        function V_tri = lam2v( this )
            % V_tri = lam2v( this )
            % Recover solution in original complete basis
            V_tri = zeros(3,3);
            for k=1:3
                V_tri(:,k) = this.Om{k} * this.lam{k};
            end
        end
        
        function v2lam( this, V_tri )
            % v2lam( this, V_tri )
            % Recover solution in original complete basis
            for k=1:3
                this.lam{k} = this.Om{k}' * V_tri(:,k);
            end
        end
        
        function V_tri = rho2v( this )
            % V_tri = rho2v( this )
            % Recover solution in original complete basis
            V_tri = zeros(3,3);
            for k=1:3
                V_tri(:,k) = this.Om{k} * this.rho(:,k);
            end
        end
        
        % Complete method
        function V_tri = solve( this )
            this.checkJunction;
%             this.getReducedJunction;
%             this.getReducedParams;
%             this.chooseTrihedronNormals;
            
            % Refine in non-diagonal basis (for most general case)
%             this.rho2lam;
%             if any(this.c ~= 0)
%                 this.lam = this.refineReducedSystem( this.B, this.c, this.lam );
%             end
%             V_tri = this.lam2v;
            
            V_tri = this.rho2v;
            this.V = V_tri;
        end
        
        function V_tri = solve_withGT( this, V_gt )
            this.getReducedParams;
%             this.chooseTrihedronNormals;
            
            % Refine in non-diagonal basis (for most general case)
            this.v2lam( V_gt );
            this.lam = this.refineReducedSystem( this.B, this.c, this.lam );
            V_tri = this.lam2v;
        end
    end
    
end