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
        
        % Auxiliar useful variables
        a0 % Product of cos among normals
        
        % Solution selection variables
        Dproj
        v_im
        Det_sign = +1 % Sign to choose for solution determinant when ambiguity is solved
        WITH_SIGNED_DIRECTIONS  = true % Flag to correct vanishing direction sign wrt segments
        WITH_SIGNED_DETERMINANT = true % Flag to choosen one solution from two possible (ambiguity)
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
            this.Dproj = repmat( {DprojectionP2( xi.c.X )}, 1,3 );
            this.v_im  = { xi.v(1).X , xi.v(2).X , xi.v(3).X };
        end
        
        function loadSegments( this, segments )
            for k=1:3
                this.Dproj{k} = DprojectionP2( segments(k).p1 );
                this.v_im{k}  = segments(k).v;
            end
            
            % Signed directions
            this.WITH_SIGNED_DIRECTIONS = true;
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
        
        function V = checkJunction( this )
            n = mat2cell(this.Nbp,3,[1 1 1]);
            
            % Check trivial case of meeting lines
            if this.isMeeting( this.Nbp )
                % Consider simple meeting lines
                IS_MEETING = true;
                
                % Other option:
%                 a  = -dot( this.Nbp, circshift(this.Nbp,[0 -1]), 1 );
%                 a0 = sqrt( a );
                a12 = -n{1}'*n{2};
                a23 = -n{2}'*n{3};
                a31 = -n{3}'*n{1};
                a0 = sqrt( a12 * a23 * a31 ); %#ok<*PROP>
                a0 = repmat( a0, 1,2 );
                a0_sign = [-1 +1];
                %                 this.getReducedJunction;
                %                 this.chooseTrihedronNormals;
                %                 return
            else
                IS_MEETING = false;
                % Compute intersection direction t (of planes 1 and 2)
                t = snormalize( cross( n{1}, n{2} ) );
                Om_t = null(t');
                
                % delta is the reduced vector for n3 direction,
                % n3 = null(t') * delta
                % Define new bilinear form (membership quadratic form)
                nm = n{3};
                if 0
                    Q = this.membership_quadratic( n{1}, n{2}, n{3} );
                else
                    M = (nm' * t)^2 / (n{1}' * skew(t)*skew(t) * n{2})...
                        * ( skew(t) * n{2} * n{1}' * skew(t) ) + nm*nm';
                    M = Om_t' * skew(t) * M * skew(t) * Om_t;
                    Q = 0.5 * (M+M'); % Symmetrize
                end
                %             disp( [Q - this.membership_quadratic( n{1}, n{2}, n{3} )] );
                
                [U,D] = eig( Q ); % Diagonalize matrix
                if sign(D(1,1)) == sign(D(2,2))
                    warning('Not different signs. Check angles condition for feasible triplets.');
                    keyboard
                    return
                end
                delta = U * snormalize(...
                    [ +sqrt(abs(D(2,2))), +sqrt(abs(D(1,1)))
                      +sqrt(abs(D(2,2))), -sqrt(abs(D(1,1))) ]');
                n3 = Om_t * delta;
                  
                % Check good solution (for each a0 sign)
                a12 = n{1}'*skew(t)*skew(t)*n{2};
                a23 = n{2}'*skew(t)*skew(t)*n3;
                a31 = n{1}'*skew(t)*skew(t)*n3;
                a0  = sqrt( [a12 a12] .* a23 .* a31 );
                err_m = -nm'*t*a0 + a12 * nm'*skew(t)*n3;
                err_p = +nm'*t*a0 + a12 * nm'*skew(t)*n3;
                err1 = [ err_m(1); err_p(1) ];
                err2 = [ err_m(2); err_p(2) ];
                % There are two a0 signs and two possible n3,
                % so 4 combinations exist but only two are feasible
                % abs(a0) is determined by n3, but a0 sign must be stored
                [~,ind(1)] = min(abs(err1));
                [~,ind(2)] = min(abs(err2));
                
                a0_sign = zeros(1,2);
                for ii=1:2
                    if ind(ii) == 1
                        a0_sign(ii) = -1;
                    elseif ind(ii) == 2
                        a0_sign(ii) = +1;
                    else
                    keyboard
                    error('Sth went wrong with sol disambiguation');
                    end
                end

%                 if norm( err_m ) < 1e-10
%                     a0_sign = [+1 +1];
%                 elseif norm( err_p ) < 1e-10
%                     a0_sign = [-1 -1];
%                 end
            end
            
%             a0_sign = [+1 -1];
            for ii = 1:2
%                 this.Nbp(:,3) = n3(:,ind(ii));
                if ~IS_MEETING
                    this.Nbp(:,3) = n3(:,ii);
                end
                
                this.getReducedJunction;
                if IS_MEETING
                    this.rho = diag([a0_sign(ii) 1]) * this.rho; % Apply 
                end
                if this.WITH_SIGNED_DIRECTIONS
                    %  Correct current rho to meet segment signs
%                     this.setSigns( diag([a0_sign(ii) 1]) * this.rho );
                    this.setSigns( diag([a0_sign(ii) 1]) * this.rho );
                end
                V_ = this.rho2v;
                % The determinant is not fixed by a0 sign, even if fixing
                % signed directions
% %                 Det = sign( det( V_ ) );
% %                 % Store det+ solution in first position and det- afterwards
% %                 if Det > 0
% %                     this.V{1} = V_;
% %                 else
% %                     this.V{2} = V_;
% %                 end
                this.V{ii} = V_;
            end
            
            if this.WITH_SIGNED_DETERMINANT
                dets = [ sign( det( this.V{1} ) ), ...
                         sign( det( this.V{2} ) ) ];
                V = this.V{ dets == this.Det_sign };
%                 if this.Det_sign > 0
%                     V = this.V{1};
%                 else
%                     V = this.V{2};
%                 end
                this.V = V;
            else
                V = this.V; % Return cell
            end
            
            % Store the final a0 value for further uses and identification
            % of solution
            this.a0 = a0 .* a0_sign;
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
                d(k) = -n{i}'*n{j};
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
            V_tri = this.checkJunction;
%             this.getReducedJunction;
%             this.getReducedParams;
%             this.chooseTrihedronNormals;
            
            % Refine in non-diagonal basis (for most general case)
%             this.rho2lam;
%             if any(this.c ~= 0)
%                 this.lam = this.refineReducedSystem( this.B, this.c, this.lam );
%             end
%             V_tri = this.lam2v;
            
%             V_tri = this.rho2v;
%             this.V = V_tri;
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
    
    methods (Static)
        % Solver functions
        function Q = membership_quadratic( n1, n2, n_ )
            t   = snormalize( cross(n1,n2) );
            Om_t = null(t');
            
            a12 = -n1'*n2;
            at_ =   t'*n_;
            
            Q = Om_t' * ...
                ( at_^2/a12 * 0.5*(n1*n2'+n2*n1') + ...
                  skew(t) * (n_*n_') * skew(t) ) * Om_t;
        end
        
        function s = signature( Nbp )
            % Signature of P3OA problem is defined as that of the
            % membership quadratic arising from three interpretation planes
            % If all three planes form a pencil the quadratic Q is
            % rank-defficient, since at_ is zero
            % If no pencil exists, Q must be indefinite for a solution to exist
            if abs( det( Nbp ) ) < 1e-10
                % Trivial case of meeting lines
                s = 0; return
            end
            Q = CTrihedronSolver.membership_quadratic( Nbp(:,1), Nbp(:,2), Nbp(:,3) );
            s = sign( det( Q ) );
        end
        
        function is = isMeeting( Nbp )
            is = ( abs( det( Nbp ) ) < 1e-12 );
        end
        
        function out = cos_prod( Nbp )
            out = prod( -dot( Nbp, circshift(Nbp,[0 -1]), 1 ) );
        end
    end
end