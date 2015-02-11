classdef CP3oASolver < handle
    % CP3oASolver Class for solving P3oA problem
    % This class has been projected as a factory of methods
    % to apply when solving tasks related to P3oA problem
    % The only property kept by class is Calibration Matrix (K)
    % since it is usually constant along a set of problems
    % The rest are just passed as input and outputs
    
    properties
        % Input
        K   % Camera intrinsic calibration
    end
    
    properties (Constant)
        % Internal properties of class
        thresDegenerate = 0.99; % Threshold below which the problem is considered degenerate
        thresDiscr = 0.01; % Threshold below which the problem can be considered degenerated (based on discriminant)
        % thresDiscr has been observed to be linked to tr(A_eps)~1 for
        % uncertainty ~1/Seg.d
        % TODO: Relate threshold to some probabilistic parameter of certainty
    end
    
    methods
        function this = CP3oASolver( K )
            if ~exist('K','var')
                error('Calibration Matrix K must be passed to CP3oASolver constructor');
            end
            this.K   = K;
        end
        
        % Auxiliar loading functions
        function loadXi( this, xi )
            this.Dproj = repmat( {DprojectionP2( xi.c.X )}, 1,3 );
            this.v_im  = { xi.v(1).X , xi.v(2).X , xi.v(3).X };
        end
                        
        % Uncertainty estimation
        function A_eps = computeCovariance( this, V, Nbp, A_Nbp )
            % Compute covariance by 1st propagation from
            % normals of interpretation planes
            if ~exist('V','var'), V = this.V; end %#ok<NASGU>
            if ~exist('Nbp','var'), Nbp = this.Nbp; end %#ok<NASGU>
            if ~exist('A_Nbp','var'), A_Nbp = this.A_Nbp; end %#ok<NASGU>
                           
            % Compute covariance of trihedron normals
            J_Phi_eps = cross( V, Nbp, 1 )';
            J_Phi_Nbp = mat2cell( V, 3, [1 1 1] );
            J_Phi_Nbp = blkdiag( J_Phi_Nbp{:} )';
            J_eps_Nbp = - J_Phi_eps \ J_Phi_Nbp;

            if numel( A_Nbp ) == 3
                % Diagonal non-correlated case
                A_N = blkdiag( A_Nbp{:} );
            elseif all( size( A_Nbp ) == [3 3] )
                % Full (correlated) case
                A_N = cell2mat( A_Nbp );
            else
                error('Wrong A_Nbp dimensions, check');
            end
            A_eps = J_eps_Nbp * A_N * J_eps_Nbp';
        end
               
        % Estimate normals covariance from segment lines covariance
        function A_N = covProp_N_L( this, A_L, segs )
            % A_N = covProp_N_L( this, segs, A_l )
            % Get normals covariance estimation from image lines covariance
            
            % Compute and collect jacobians for each line
            % Different lines are independent of each other
            Jac_n_l = cell(1,3);
            for k=1:3
                % Compute jacobian for each line
                Jac_n_l{k} = Dsnormalize(this.K'*segs(k).l) * this.K';
            end
            
            if numel( A_L ) == 3
                % Diagonal non-correlated case (lines are independent)
                A_N = cell(1,3);
                for k=1:3
                    A_N{k} = Jac_n_l{k} * A_L{k} * Jac_n_l{k}';
                end
            elseif all( size( A_L ) == [3 3] )
                % Full (correlated) case
                Jac_N_L = blkdiag( Jac_n_l{:} );
                A_L     = cell2mat( A_L );
                A_N     = Jac_N_L * A_L * Jac_N_L';
                A_N = mat2cell(A_N,[3 3 3],[3 3 3]); % Divide in cell again
            else
                warning('Wrong A_L dimensions, check');
                keyboard
            end
        end
        
        % Use uncertainty propagation to create ficticial weight
        function A_eps = covProp_eps_N( this, A_N, Nbp, V )
            % Compute covariance of triplet directions
            J_Phi_eps = cross( V, Nbp, 1 )';
            J_Phi_Nbp = mat2cell( V, 3, [1 1 1] );
            J_Phi_Nbp = blkdiag( J_Phi_Nbp{:} )';
            J_eps_Nbp = - J_Phi_eps \ J_Phi_Nbp;
            
            if numel( A_N ) == 3
                % Diagonal non-correlated case (normals are independent)
                A_N = blkdiag( A_N{:} );
            elseif all( size( A_N ) == [3 3] )
                % Full (correlated) case
                A_N = cell2mat( A_N );
            else
                warning('Wrong A_N dimensions, check');
                keyboard
            end
            
            % Propagation step
            A_eps = J_eps_Nbp * A_N * J_eps_Nbp';
        end
        
        function V = solveMeetingCase( this, Nbp, necker_sign )
            % V = solveMeetingCase( Nbp, necker_sign )
            % Compute solution for the meeting lines case
            % Necker sign stands for solution duality
            % according to Necker cube illusion
            
            n = mat2cell(Nbp,3,[1 1 1]); % Interpretation normals
            % Compute interpretation cosines
            % Operation below stands for:
            a12 = -n{1}'*n{2};
            a23 = -n{2}'*n{3};
            a31 = -n{3}'*n{1};
%             a  = -dot( this.Nbp, circshift(Nbp,[0 -1]), 1 );
            % a0 is the abs value of a_bar parameter, which must
            % be complete with Necker sign
            %   a_bar = +- sqrt( a12*a23*a31 )
%             a_bar = necker_sign * sqrt( prod( a ) );
            a_bar = necker_sign * sqrt( a12*a23*a31 );
            if ~isreal(a_bar)
                warning('Complex solution, check angles condition');
                keyboard
            end
            
            % Apply closed solution to reduced lambda parameter
            % Use specified Necker sign
%             lam = snormalize([ a_bar a_bar a_bar ; a ]);
            lam = snormalize([ a_bar a_bar a_bar ;
                                a23   a31   a12  ]);
            
            % Define conversion bases (lam to v)
            t = snormalize( cross( n{1}, n{2} ) ); % Vertex ray
            Om{1} = [ t, cross(t,n{1}) ];
            Om{2} = [ t, cross(t,n{2}) ];
            Om{3} = [ t, cross(t,n{3}) ];
            
            V = zeros(3,3);
            for k=1:3
                V(:,k) = Om{k} * lam(:,k);
            end
        end

        % Complete method
        function cell_V = solve( this, Nbp )
            % cell_V = solve( Nbp )
            % Solve the general problem of P3oA
            % given three interpretation planes Nbp (three normals)
            
            cell_V = cell(2,1);
            if this.isMeeting( Nbp )
                % Simple meeting lines case
                cell_V{1} = this.solveMeetingCase( Nbp, +1 );
                cell_V{2} = this.solveMeetingCase( Nbp, -1 );
            else
                % General no meeting lines case
                n = mat2cell(Nbp,3,[1 1 1]);
                t = snormalize( cross( n{1}, n{2} ) );
                Om_t = null(t');
                
                % Compute membership quadratic form
%                 Q = this.membership_quadratic( n{1}, n{2}, n{3} );
                nm = n{3};
                M = (nm' * t)^2 / (n{1}' * skew(t)*skew(t) * n{2})...
                    * ( skew(t) * n{2} * n{1}' * skew(t) ) + nm*nm';
                M = Om_t' * skew(t) * M * skew(t) * Om_t;
                Q = 0.5 * (M+M'); % Symmetrize
                
                % Solve the reduced quadratic problem as eigenvalue problem
                [U,D] = eig( Q ); % Diagonalize matrix
                if sign(D(1,1)) == sign(D(2,2))
                    warning('Wrong Q. Signature should be [0 1 1]');
                    keyboard
                    return
                end
                % In the diagonal representation one can see two signs
                % are undetermined so 4 possible solutions exist for delta
                % However, only two are linearly independent
                % We take 
                delta = zeros(2,2);
                sqrtD = sqrt(abs(diag(D)));
                delta(:,1) = U * snormalize( [ +sqrtD(2); +sqrtD(1) ] );
                delta(:,2) = U * snormalize( [ +sqrtD(2); -sqrtD(1) ] );
                % Sign is changed in the biggest element D(1,1) to have
                % greater difference between solutions
                
                % Compute two virtual n3 vectors which complete the meeting case
                n3_ = Om_t * delta;
                
                % Define 4 possible solutions (cartesian product of
                % Necker sign and virtual n3_)
                sol = repmat( struct('nvir',[],'necker_sign',[],'err',[]), 4,1 );
                sol(1).nvir = n3_(:,1);
                sol(1).necker_sign = +1;
                sol(2).nvir = n3_(:,2);
                sol(2).necker_sign = +1;
                sol(3).nvir = n3_(:,1);
                sol(3).necker_sign = -1;
                sol(4).nvir = n3_(:,2);
                sol(4).necker_sign = -1;
                % Check membership condition error for each solution
                for k=1:4
                    nvir = sol(k).nvir;
                    necker_sign = sol(k).necker_sign;
                                        
                    a12 = -n{1}'*n{2};
                    a23 = -n{2}'*nvir;
                    a31 = -nvir'*n{1};
                    a_bar = necker_sign * sqrt( a12*a23*a31 );
                    if ~isreal(a_bar)
                        warning('Complex solution, check angles condition');
                        keyboard
                    end
                    t = snormalize( cross( n{1}, n{2} ) ); % Vertex ray
                    err = a_bar * n{3}'*t + a12 * n{3}'*skew(t)*nvir;
                    
                    % Store membership error
                    sol(k).err = abs(err);
                end
                % Sort errors to find the two good solutions
                [~,I] = sort([sol.err],2,'ascend');
                
                % Apply simple meeting lines case to good solutions
                for k=1:2
                    Nbp = [n{1},n{2},sol(I(k)).nvir];
                    necker_sign = sol(I(k)).necker_sign;
                    cell_V{k} =...
                        this.solveMeetingCase( Nbp, necker_sign );
                end
            end
        end
        
        % Method to correct direction sign wrt given oriented segments
        function V = correctSigns( this, V, segs )
            % V = correctSigns( this, segs )
            % Receives a vector of oriented segments
            % corresponding to the three directions in V
            % Return corrected V so that its vectors project
            % with proper sign
            
            for k=1:3
                v_im  = segs(k).v; % Oriented vector in image
                Dproj = DprojectionP2( segs(k).p1 ); % Proj derivative in p1 point
                
                % Correct V k-th direction with dot product sign
                s = sign( v_im' * Dproj * this.K * V(:,k) );
                V(:,k) = s * V(:,k);
            end
        end
        
        % Method to choose good solution wrt given oriented segments and
        % determinant sign
        function V = correctSignsDet( this, cell_V, segs, detSign )
            % V = correctSignsDet( this, segs, detSign )
            % Receives:
            % - a vector of oriented segments
            % corresponding to the three directions in V
            % - a sign for determinant of set of directions
            % Return corrected V so that its vectors project
            % with proper sign and the system has given determinant sign
            %
            % Note: after applying both conditions, only one solution
            % stands
            
            for k=1:2
                V = this.correctSigns( cell_V{k}, segs );
                if sign(det(V))==detSign
                    return
                end
            end
            
            % Check complex values
            if ~isreal( V )
                warning('Complex values, not valid problem');
                V = NaN(3);
                return
            end
            
            % If not returned from inside for-loop, sth went wrong
            warning('No detSign found? Check');
%             keyboard
            V = NaN(3);
        end
    end
    
    methods (Static)
        % Solver functions
        % General method unifying all cases
        function cell_V = generalSolve( Nbp )
            % Alternative solution using eigendecomposition
            
            n = mat2cell(Nbp,3,[1 1 1]);
            
            % Compute interpretation planes bases
            % NOTE: Pencil bases can be used to
            % simplify the numeric values as much as possible
            Om = cell(1,3);
            for k=1:3
                Om{k} = null(n{k}');
            end
            
            % Compute bilinear forms
            B = cell(1,3);
            for cidx=1:3
                circ_idxs = num2cell( circshift(1:3,[0 cidx]) );
                [i,j,k] = deal(circ_idxs{:});
                B{k} = Om{i}'*Om{j};
            end
            % 90 degree 2D rotation or orthogonal rotation
            T = [0 -1
                 1  0];
            
            % Solve eigendecomposition of combined matrices          
            cell_V = {zeros(3);zeros(3)};
            for cidx=1:3
                circ_idxs = num2cell( circshift(1:3,[0 cidx]) );
                [i,j,k] = deal(circ_idxs{:});
                % Compute concatenation matrix
                M = T*B{k}*T*B{i}*T*B{j};
                [P,D] = eig( M );
                % Sort eigenvalues so that solutions to two Necker
                % configurations are properly collected
                % Eigenvalues sometimes are not returned ordered
                [~,I] = sort( diag(D) );
                cell_V{1}(:,i) = Om{i} * P(:,I(1));
                cell_V{2}(:,i) = Om{i} * P(:,I(2));
            end
        end
        
        function A = discriminant( Nbp )
            % A = discriminant( Nbp )
            % In the general solution, compute discriminant of
            % characteristic polynomial (2nd order) of any
            % concatenation matrix to check if real solutions exist
            
            if 1
                % Compute interpretation cosines
                % IMPORTANT: sign is contrary of interpretation cosine a
                % which are aij = -ni'*nj
                c = dot( Nbp, circshift(Nbp,[0 -1]), 1 );
                
                % Compute det of concatenation matrix
                % Since M is 2x2 and its determinant and trace
                % are invariant to cyclic permutation, 
                % is characteristic polynomial
                % (which only depends on tr(M) and det(M))
                % is the same for all 3 cases of concatenation matrix M
                % det(M)   = (n1'*n2) * (n2'*n3) * (n3'*n1)
                % trace(M) = det([n1 n2 n3])
                DetM = prod(c);
                TrM  = det(Nbp);
                A = TrM^2 - 4*DetM;
                
            else
                % DEPRECATED:
                % Other option is algebraically optimized
                % to avoid all non-necessary operations! :D
                % Extract normals
                n = mat2cell(Nbp,3,[1 1 1]);

                % TODO: Reorder code to avoid repeating so many operations

                % Compute interpretation planes bases
                % NOTE: Pencil bases can be used to
                % simplify the numeric values as much as possible
                Om = cell(1,3);
                for k=1:3
                    Om{k} = null(n{k}');
                end

                % Compute bilinear forms
                B = cell(1,3);
                for cidx=1:3
                    circ_idxs = num2cell( circshift(1:3,[0 cidx]) );
                    [i,j,k] = deal(circ_idxs{:});
                    B{k} = Om{i}'*Om{j};
                end
                % 90 degree 2D rotation or orthogonal rotation
                T = [0 -1
                     1  0];

                M = T*B{1}*T*B{2}*T*B{3};
                p = charpoly(M); % Characteristic polynomial
                A = p(2)^2 - 4*p(1)*p(3);
            end
        end
        
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
            Q = CP3oASolver.membership_quadratic( Nbp(:,1), Nbp(:,2), Nbp(:,3) );
            s = sign( det( Q ) );
        end
        
        function is = isMeeting( Nbp )
            is = ( abs( det( Nbp ) ) < 1e-12 );
        end
        
        function out = cos_prod( Nbp )
            out = prod( -dot( Nbp, circshift(Nbp,[0 -1]), 1 ) );
        end
        
        function [is, d] = isDegenerate( Nbp )
            % [is, d] = isDegenerate( Nbp )
            % Return true if the configuration is degenerate
            % An internal threshold (thresDegenerate) will be used
            % to decide if the problem is degenerate based on the
            % interpretation cosines 

            if 0
                % Compute interpretation cosines
                c = -dot( Nbp, circshift(Nbp,[0 -1]), 1 );

                % Check that no pair of normals is above maximum dot similarity
                d = max( abs(c) );
                is = (d > CP3oASolver.thresDegenerate);
            else
                d  = CP3oASolver.discriminant( Nbp );
                is = ( d < CP3oASolver.thresDiscr );
            end
        end
    end
end
