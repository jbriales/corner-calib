classdef CTrihedronSolver < handle
    %CTrihedronSolver Class for solving trihedron problem
    %   Detailed explanation goes here
    
    properties
        % Input
        Nbp % Back-projected planes normals
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
        
        % Solver functions
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
            
% % %             % Linear system
% % %             M = cell(6,6);
% % %             M{1,5} = U{1}; M{2,6} = U{1};
% % %             M{1,3} = -V{1}; M{2,4} = -V{1};
% % %             M{3,1} = U{2}; M{4,2} = U{2};
% % %             M{3,5} = -V{2}; M{4,6} = -V{2};
% % %             M{5,3} = U{3}; M{6,4} = U{3};
% % %             M{5,1} = -V{3}; M{6,2} = -V{3};
% % %             for ii=1:numel(M)
% % %                 if isempty(M{ii})
% % %                     M{ii} = zeros(2,2);
% % %                 end
% % %             end
% % %             MM = cell2mat( M );
% % %             NN = null(MM);
% % %             % Temporal
% % %             WWa = mat2cell(reshape(NN(:,1),2,[]), 2, [2 2 2]);
% % %             WWb = mat2cell(reshape(NN(:,2),2,[]), 2, [2 2 2]);
% % %             [Wa{1},Wa{2},Wa{3}] = deal( WWa{:} );
% % %             [Wb{1},Wb{2},Wb{3}] = deal( WWb{:} );
% % %             
% % %             % Check random transformations
% % %             for k=1:3
% % %                 Q{k} = rand(2);
% % %             end
% % %             [U_{1},S_{3},V_{2}] = svd( Q{1}'*B{1,2}*Q{2} );
% % %             [U_{2},S_{1},V_{3}] = svd( Q{2}'*B{2,3}*Q{3} );
% % %             [U_{3},S_{2},V_{1}] = svd( Q{3}'*B{3,1}*Q{1} );
% % %             
% % %             [vv,dd] = eig( B{1,2} )
% % %             [vv,dd] = eig( B{2,3} )
% % %             [vv,dd] = eig( B{3,1} )
% % %             
% % %             % Seek simultaneous singular eigenvectors
% % %             A{3} = inv(B{1,3}) * inv(B{3,1}) * ...
% % %                    B{3,2} * B{2,3};
% % %             A{3} = B{3,1} * B{1,3} * B{3,2} * B{2,3};
% % %             A{2} = B{2,1} * B{1,2} * B{2,3} * B{3,2};
% % %             A{1} = B{1,2} * B{2,1} * B{1,3} * B{3,1};
% % %             [Q{3},~,T{3}] = svd( A{3} );
% % %             [Q{2},~,T{2}] = svd( A{2} );
% % %             [Q{1},~,T{1}] = svd( A{1} );
% % %             
% % %             for k=1:3
% % %                 Om{k} = Om{k} * Q{k};
% % %             end
% % %             for k=1:3
% % %                 ij = setdiff(1:3,k); i=ij(1); j=ij(2);
% % %                 B{i,j} = Om{i}' * Om{j};
% % %                 B{j,i} = Om{j}' * Om{i};
% % %             end
% % %             
% % %             % Big bilinears
% % %             Os = zeros(2);
% % %             BB{3} = 0.5 * [ Os B{1,2} Os
% % %                             B{1,2} Os Os
% % %                             Os Os Os ];
% % %             BB{2} = 0.5 * [ Os Os B{1,3}
% % %                             Os Os Os
% % %                             B{1,3} Os Os ];
% % %             BB{1} = 0.5 * [ Os Os Os
% % %                             Os Os B{2,3}
% % %                             Os B{2,3} Os ];
% % %             
% % %             % Diagonalize bilinear forms
% % %             U = cell(1,3);
% % %             [U{1},S{3},V{2}] = svd( B{1,2} );
% % %             [U{2},S{1},V{3}] = svd( B{2,3} );
% % %             [U{3},S{2},V{1}] = svd( B{3,1} );
% % %             
% % %             % Create new system for nullspace vector computation
% % %             syms a b
% % %             for k=1:3
% % %                 W_{k} = a*Wa{k} + b*Wb{k};
% % %             end
% % %             for k=1:3
% % %                 eq{k} = W_{k}(:,1).' * S{k} * W_{k}(:,2);
% % %             end
            
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
                              this.K * this.OmU{i} *...
                              rho(:,i) );
                    rho(:,i) = s * rho(:,i);
                    n{i} = this.OmU{i} * rho(:,i);
                end
                R_tri = cell2mat( n );
                
                if abs( det(R_tri) - 1 ) < 0.5
                    break % Exit loop with current rho output value
                end
            end
            
            this.rho = rho; % Store good solution
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
                V_tri(:,k) = this.OmU{k} * this.rho(:,k);
            end
        end
        
        % Complete method
        function V_tri = solve( this )
            this.getReducedParams;
            this.chooseTrihedronNormals;
            
            % Refine in non-diagonal basis (for most general case)
            this.rho2lam;
            this.lam = this.refineReducedSystem( this.B, this.c, this.lam );
            V_tri = this.lam2v;
            
%             V_tri = this.rho2v;
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