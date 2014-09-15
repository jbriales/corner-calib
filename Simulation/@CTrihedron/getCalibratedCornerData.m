function [N, c, A_co, L_P2, A_L_P2] = getCalibratedCornerData( obj, img_pts, Camera )
% [N, c, A_co] = getCalibratedCornerData( obj, Camera )
% Input:
%   Camera  - CSimCamera object grabbing the pattern
%
% Output:
%   img     - img structure with the additional fields
%       x   - 5xN array with: c point (in pixels), 3 angles of direction
%       A_x - 5x5xN covariance matrices with the uncertainty of x

% Control variables
sigma = Camera.sd; % TODO: Square? sigma = sd^2?

% Define the x and A_x structures
x   = zeros(5,1);
A_x = zeros(5,5,1);

% Assignment of auxiliar variables
c   = img_pts(:,1);
p   = img_pts(:,2:4) - repmat(c,1,3);
a   = atan2( p(2,:), p(1,:) );

% Estimation of x (unused)
x(1:2) = c;
x(3:5) = a;
% Estimation of A_x
A_x(1,1) = sigma;
A_x(2,2) = sigma;

A_pc = sigma * eye(4);
J    = zeros(3,4);
for i = 1:3
    % Set the num and den variables
    px = p(1,i);
    py = p(2,i);
    % Set the derivatives of atan2
    alpha = px.^2 + py.^2;
    d_atan_x = -py./alpha;
    d_atan_y =  px./alpha;
    % Calculate the Jacobian matrix
    J(i,1) =  d_atan_x;
    J(i,2) =  d_atan_y;
    J(i,3) = -d_atan_x;
    J(i,4) = -d_atan_y;
end
A_x(3:5,3:5) = J * A_pc * J';

%% Convert to calibrated data
K = Camera.K;
% Compute homogeneous lines
l = cell(1,3);
for k=1:3
    n = [-sin(a(k)) cos(a(k))]';
    d = - n' * c;
    l{k} = [ n ; d ];
end

% Calibrate image
% Central point and directions are converted from pixel values to m values
J_an_a = zeros(3,3);
for k=1:3
    n    = l{k}(1:2); % Normal in pixels space
    
    l{k} = K' * l{k};
    l{k} = l{k} / norm(l{k}(1:2));
    
    g = K(1:2,1:2)' * n; % Auxiliar variable: non-normalized normal in [m]
    n_n = g / norm(g); % Normal in [m] space
    
    % Compute the Jacobians of the normal directions
    J_an_nn = [-n_n(2) n_n(1)];
    J_nn_g  = [g(2)^2 -g(1)*g(2);-g(1)*g(2) g(1)^2]/norm(g)^3;
    J_g_n   = K(1:2,1:2)';
    J_n_a   = [-n(2) n(1)]';
    J_an_a(k,k) = J_an_nn * J_nn_g * J_g_n * J_n_a;
end
% Compute the Jacobians of the propagated point
K_inv   = inv(K);
c_n     = K \ makehomogeneous(c);
J_cn_c  = [1/c_n(3) 0 -c_n(1)/c_n(3)^2;
    0 1/c_n(3) -c_n(2)/c_n(3)^2 ] * K_inv(1:3,1:2);
c = makeinhomogeneous( c_n );

%% Output data
% Set uncertainty:
J    = [J_cn_c zeros(2,3); zeros(3,2) J_an_a];
% J    = [J_an_a zeros(3,2); zeros(2,3) J_cn_c];
A_co = J * A_x * J';

% Calibrated homogeneous lines
L_P2 = cell2mat( l );
J_L_CO = Jacobian_L_CO( c, l{1}(1:2), l{2}(1:2), l{3}(1:2) );
A_L_P2 = J_L_CO * A_co * J_L_CO';

% TODO: Compute homogeneous lines uncertainty (and also for reprojection
% planes normals)

% Minimal data for rotation
N = L_P2(1:2,:);

end

function [ J_L_CO ] = Jacobian_L_CO ( c, n1, n2, n3 )

Jc = @(n) [ zeros(2,2) ; -n' ];
Jn = [ eye(2) ; -c' ];
J_L_cN = [ [Jc(n1); Jc(n2); Jc(n3)] , blkdiag(Jn,Jn,Jn) ];

T = [0 -1 ; 1 0];
J_cN_co = blkdiag( eye(2), T*n1, T*n2, T*n3 );

J_L_CO = J_L_cN * J_cN_co;

end