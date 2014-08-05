function [N, A_N, A_eps] = getWorldNormals( R0, L, c, A_CO )
% [N, AN] = getWorldNormals( L, c, AL, Ac )
% Function that receives Corner Observation (and its covariance matrix) and
% computes the world normals (and its covariance matrix too)
% Input:
%   R0 is the initial estimate for rotation matrix R_c_w
%   L is a 2x3 array of line normals (in 2D)
%   c is a 2x1 column vector with corner point
%   A is the covariance matrix (5x5) of normals (1D) and corner point (2D)

% % Extract covariance components
% A_l1 = A_CO(1,1);
% A_l2 = A_CO(2,2);
% A_l3 = A_CO(3,3);
% A_c  = A_CO(4:5,4:5);

% Extract Corner data
l1 = L(:,1);
l2 = L(:,2);
l3 = L(:,3);

%% Solve rotation estimation
% Coefficient matrix
M = [ l1' -l1'*c
      l2' -l2'*c
      l3' -l3'*c ];

% Implicit function to solve in SO(3) space and its Jacobian for Newton
% optimization
% R is R_c_w (World orientation seen from Camera)
% N = [nx ny nz] is a 3x3 array with normal directions tiled as column
% vectors
F_Phi = @(R) dot( M', R, 1 )';
J_Phi = @(R)[ -M(1,:)*skew(R(:,1))
              -M(2,:)*skew(R(:,2))
              -M(3,:)*skew(R(:,3)) ];
          
debug = 0;
R = Newton.optim(F_Phi, J_Phi, R0,'SO(3)', debug);

N = R;

% Uncertainty propagation
% Check dimensions
if size(A_CO,1) ~= 5
    error('A_CO dimensions should be 5 (1+1+1+2)')
end

% Debug
% A_CO = 1e-4 * eye(5); % TO REMOVE

[ J_R_CO, J_eps_CO ] = Jacobian_R_CO(c,l1,l2,l3,R);
A_R = J_R_CO * A_CO * J_R_CO';
A_eps = J_eps_CO * A_CO * J_eps_CO';
A_N = A_R;

global WITH_MONTECARLO
if WITH_MONTECARLO
%     MonteCarlo_simulate( c, l1, l2, l3, A_CO, A_R, R )
    MonteCarlo_simulate( c, l1, l2, l3, A_CO, A_eps, R )
end

end

function [ J_R_CO, J_eps_CO ] = Jacobian_R_CO ( c, n1, n2, n3, R )
% Returns Jacobian of rotation (tangent space) wrt to data in corner
% observation

% Jacobian of coefficients matrix M wrt directions (1, 2 and 3): 3x2
J_M_n = [ eye(2)
             -c' ];
         
% Set Jacobian in S1 tangent space (for line directions)
% Jacobian of coefficients matrix M wrt corner point: 9x2
J_M_c = [ zeros(2)
              -n1'
          zeros(2)
              -n2'
          zeros(2)
              -n3' ];
% IMP: Use this step only in input covariance matrix is scalar for
% directions (minimal representation)
J_M_n1 = J_M_n * [ 0 -1 ; 1 0 ] * n1;
J_M_n2 = J_M_n * [ 0 -1 ; 1 0 ] * n2;
J_M_n3 = J_M_n * [ 0 -1 ; 1 0 ] * n3;

% Complete Jacobian of M wrt Corner Observation (n1,n2,n3,c)
J_M_CO = [ J_M_c, blkdiag( J_M_n1, J_M_n2, J_M_n3 ) ];

% Implicit function theorem: Jacobian of rotation in tangent space so(3)
% wrt M
M = [ n1' -n1'*c
      n2' -n2'*c
      n3' -n3'*c ];
J_Phi_eps = [ -M(1,:)*skew(R(:,1))
              -M(2,:)*skew(R(:,2))
              -M(3,:)*skew(R(:,3)) ];

J_Phi_M   = blkdiag( R(:,1)', R(:,2)', R(:,3)' );

J_eps_M   = - J_Phi_eps \ J_Phi_M;

% Jacobian of complete SO(3) space wrt tangent space so(3)
J_R_eps  = [ -skew(R(:,1))
             -skew(R(:,2))
             -skew(R(:,3)) ];
J_eps_R  = 0.5 * J_R_eps';

% Chain rule: Jacobian of rotation in tangent space so(3) wrt Corner
% Observation
J_eps_CO = J_eps_M * J_M_CO;
J_R_CO = J_R_eps * J_eps_M * J_M_CO;

% Out:
% J = J_R_CO;

end

function MonteCarlo_simulate( c, l1, l2, l3, A_CO, A_R, R0 ) 

%% From oa,ob,oc and o to R
X0 = [c',l1',l2',l3']';

% Formation of M data from corner observation
% M = [ n1' , -n1'*o
%       n2' , -n2'*o
%       n3' , -n3'*o ]';

% F_M - function that returns M matrix from corner observation vector
% F_M = @(X) [ X(3:4,:)' , -X(3:4,:)'*X(1:2,:)
%              X(5:6,:)' , -X(5:6,:)'*X(1:2,:)
%              X(7:8,:)' , -X(7:8,:)'*X(1:2,:) ];
F_M = @(X) [ X(3:4)' , -X(3:4)'*X(1:2)
             X(5:6)' , -X(5:6)'*X(1:2)
             X(7:8)' , -X(7:8)'*X(1:2) ];

% Function Phi and its Jacobian depending on R and M
F_Phi = @(R,M) dot( M', R, 1 )';
J_Phi = @(R,M)[ -M(1,:)*skew(R(:,1))
                -M(2,:)*skew(R(:,2))
                -M(3,:)*skew(R(:,3)) ];

% F_R - function that from set of data (corner observation) returns
% rotation optimization by Newton method
F_R = @(X0) Newton.optim( @(R)F_Phi(R,F_M(X0)), @(R)J_Phi(R,F_M(X0)), R0, 'SO(3)' );

fprintf('Montecarlo analysis of R (3x3) solving from Corner Observation\n')
% tic, [mu,A] = MonteCarlo.simulate( F_R, X0, 'propagation', true, 'AF', A_R, 'A0', A_CO, 'N', 1e3, 'inSpace', {'S1','S1','S1','Rn'},'outSpace', 'SO(3)' ); toc
[manSum,manSubtraction,manMean] = MonteCarlo.spaceOps( {'R','R','S1','S1','S1'}, {'SO(3)'} );
tic, [mu,A] = MonteCarlo.simulate( F_R, X0, A_CO, 'manSum',manSum, 'manSubtraction',manSubtraction, 'manMean',manMean,...
                                   'AF', A_R, 'N', 1e3 ); toc
keyboard

% function X = manSum(X0,inc_eps)
%     NX = length(X0);
%     N  = size(inc_eps,2);
%     
%     X = zeros( NX, N );
%     X(1:2,:) = MonteCarlo.manSumRn( X0(1:2),inc_eps(1:2,:) );
%     X(3:4,:) = MonteCarlo.manSumS1( X0(3:4),inc_eps(3,:) );
%     X(5:6,:) = MonteCarlo.manSumS1( X0(5:6),inc_eps(4,:) );
%     X(7:8,:) = MonteCarlo.manSumS1( X0(7:8),inc_eps(5,:) );
% %     funs = {@MonteCarlo.manSumRn,...
% %             @MonteCarlo.manSumS1,...
% %             @MonteCarlo.manSumS1,...
% %             @MonteCarlo.manSumS1};
% % 	spans_spa = {1:2,3:4,5:6,7:8};
% %     spans_man = {1:2,3,4,5};
% %     for i=1:length(funs)
% %         idx_spa = spans_spa{i};
% %         idx_man = spans_man{i};
% %         fun = funs{i};
% %         X(idx_spa,:) = fun( X0(idx_spa), inc_eps(idx_man,:) );
% %     end
% end

end

