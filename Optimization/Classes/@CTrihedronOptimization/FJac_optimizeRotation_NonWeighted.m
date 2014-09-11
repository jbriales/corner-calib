function J = FJac_optimizeRotation_NonWeighted( obj, R )
% Compute Jacobian of R wrt input parameters n and l by implicit function
% theorem
% The size will be 3x(5N) where N is the number of correspondences
n_cam = obj.cam_N;
v_LRF = obj.LRF_V;

N = size(n_cam,2); % Nr of pairs
% H = zeros(3,3);
% D = zeros(3,5*N);
blockH = cell(1,N);
blockD = cell(1,N);
% Compute Hessian and 2nd der
for i=1:N % For each pair
    ni  = n_cam(:,i); % Normal vector seen from camera (3x1)
    lsi = v_LRF(:,i); % Line vector seen from scanner (2x1)
    lci = R(:,1:2) * lsi; % Line vector seen from camera (3x1)
    
    % Phi is dC/deps = -Sum( cross(R*lsi, ni*(R*lsi)'*ni) )
    d_phi_ni  = 2 * skew(lci) * (ni*lci' + lci'*ni*eye(3));
    d_phi_lci = 2 * skew(lci) * (ni*ni') - 2 * skew(ni * lci' * ni);
    d_lci_lsi = R(:,1:2);
    d_lci_eps = - skew( R(:,1:2)*lsi );
    
    % Derivation wrt l Lie algebra instead of complete space
    d_lsi_alpha = [ -lsi(2), +lsi(1) ]';
    
%     H = H + d_phi_lci * d_phi_eps;
    blockH{i} = d_phi_lci * d_lci_eps;
%     blockD{i} = [ d_phi_ni , d_phi_lci * d_lci_lsi ];
    blockD{i} = [ d_phi_ni , d_phi_lci * d_lci_lsi * d_lsi_alpha ];
end
H = sum( reshape(cell2mat( blockH ),3,3,[]), 3 );
D = cell2mat( blockD );

J = - H \ D; % Jacobian of rotation in Lie algebra wrt all inputs ni and li
end