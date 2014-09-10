function [N, A_N, A_eps] = getWorldNormals( obj, R0, N, c, A_CO )
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
n1 = N(:,1);
n2 = N(:,2);
n3 = N(:,3);

%% Solve rotation estimation
% Coefficient matrix
M = [ n1' -n1'*c
      n2' -n2'*c
      n3' -n3'*c ];

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

[ J_R_CO, J_eps_CO ] = Jacobian_R_CO(c,n1,n2,n3,R);
A_R = J_R_CO * A_CO * J_R_CO';
A_eps = J_eps_CO * A_CO * J_eps_CO';
A_N = A_R;

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

