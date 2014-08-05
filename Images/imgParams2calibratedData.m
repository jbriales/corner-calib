function [N, p, A_co, L_P2] = imgParams2calibratedData( x, K, A_x )
% [N, p, A_co, hL] = pts2calibratedData( P, K )
% Input:
%   x - 1x5 array with: c point (in pixels), 3 angles of direction
%   K - intrinsic calibration matrix 
% Output:
%   N  - 2x3 array with direction of lines
%   p  - 2D corner point
%   A_co - covariance matrix with uncertainty of output data
%   L_P2 - 3x3 array whose columns are (hom.) calibrated lines

%% Set input values
p = x(1:2);
a(1) = x(3);
a(2) = x(4);
a(3) = x(5);

%% Compute homogeneous lines
l = cell(1,3);
for k=1:3
    n = [-sin(a(k)) cos(a(k))]';
    d = - n' * p;
    l{k} = [ n ; d ];
end

%% Calibrate image
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
c_n     = K \ makehomogeneous(p);
J_cn_c  = [1/c_n(3) 0 -c_n(1)/c_n(3)^2;
           0 1/c_n(3) -c_n(2)/c_n(3)^2 ] * K_inv(1:3,1:2);
p = makeinhomogeneous( c_n );

%% Output data
% Calibrated homogeneous lines
L_P2 = cell2mat( l );

% Minimal data for rotation
N = L_P2(1:2,:);

% Set uncertainty:
J    = [J_cn_c zeros(2,3); zeros(3,2) J_an_a];
% J    = [J_an_a zeros(3,2); zeros(2,3) J_cn_c];
A_co = J * A_x * J';

global WITH_MONTECARLO
if WITH_MONTECARLO
    Fun = @(X)F(X(1:2),X(3:4),X(5:6),X(7:8), K);
    X0  = [x(1:2,1) ; [-sin(a(1)) cos(a(1))]' ; [-sin(a(2)) cos(a(2))]' ; [-sin(a(3)) cos(a(3))]'];
    A0  = A_x;
    [manSum,manSubtraction,manMean] = MonteCarlo.spaceOps( {'R','R','S1','S1','S1'}, {'R','R','S1','S1','S1'} );
    MonteCarlo.simulate( Fun, X0, A0, 'manSum',manSum, 'manSubtraction',manSubtraction, 'manMean',manMean,...
                                  'AF', A_co, 'N', 1e4 );
    keyboard
end

end

function out = F(p, n1, n2, n3, K)
K12 = K(1:2,1:2);
n1 = K12' * n1; n1 = n1/norm(n1);
n2 = K12' * n2; n2 = n2/norm(n2);
n3 = K12' * n3; n3 = n3/norm(n3);
p  = makeinhomogeneous( K \ makehomogeneous( p ) );
out = [ p ; n1 ; n2 ; n3 ];
end