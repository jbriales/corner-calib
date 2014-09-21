function obj_Rtri = getTrihedronNormals( obj, obj_Nbp, R0 )
% obj_Rtri = getTrihedronNormals( obj, obj_Nbp, R0 )
% Function that receives Corner Observation (and its covariance matrix) and
% computes the world normals (and its covariance matrix too)
% Input:
%   R0 is the initial estimate for rotation matrix R_c_w

% Implicit function to solve in SO(3) space and its Jacobian for Newton
% optimization
% R is Rtri (trihedron orientation seen from Camera)
% N = [nx ny nz] is a 3x3 array with normal directions tiled as column
% vectors
N = obj_Nbp.arr;
F_Phi = @(R) dot( N, R, 1 )';
J_Phi = @(R) cross( R, N, 1 )';
          
debug = 0;
Rtri = Newton.optim(F_Phi, J_Phi, R0,'SO(3)', debug);

obj_Rtri = Manifold.SO3( Rtri );

% Compute covariance of trihedron normals
J_Phi_eps = J_Phi( Rtri );

J_Phi_Nbp = mat2cell( Rtri, 3, [1 1 1] );
J_Phi_Nbp = blkdiag( J_Phi_Nbp{:} )';

J_eps_Nbp = - J_Phi_eps \ J_Phi_Nbp;
% J_Rtri_eps = obj_Rtri.Dexp;
% J_Rtri_Nbp = J_Rtri_eps * ( - J_phi_eps \ J_Phi_Nbp );

obj_Rtri.setMinimalCov( J_eps_Nbp * obj_Nbp.A_X * J_eps_Nbp' );
end