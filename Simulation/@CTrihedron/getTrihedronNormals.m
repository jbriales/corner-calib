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
% Rgt = R0;
% figure, hold on, rotate3d on, axis equal
% plotframe(eye(4), 1, 'W', 'k');
% plotframe([N zeros(3,1); zeros(1,3) 1], 0.5, 'N', 'k', 'LineWidth',3);
% plotframe([Rgt zeros(3,1); zeros(1,3) 1], 0.5, 'gt', 'rgb', 'LineWidth',3);
% for i=1:50
%     [U,~,V] = svd( randn(3,3) );
%     R0 = U*V' * diag([det(U*V') 1 1]);
% %     R0 = ( Rgt' * RotationZ( deg2rad( 30*rand(1) ) ) *...
% %            RotationX( 0.3*rand(1) ) * RotationY( 0.3*rand(1)) )';
% 
%     Rtri = Newton.optim(F_Phi, J_Phi, R0,'SO(3)', debug);
%     fprintf('Distance: %f\n',angularDistance(Rgt,Rtri));
%     disp([R0 Rgt Rtri])
%     h = plotframe([Rtri zeros(3,1); zeros(1,3) 1], 0.5);
%     keyboard
%     delete(h)
% end

% Plot Phi function
if 0
    F_Cost = @(R) squeeze( sum( dot( repmat(N,[1 1 size(R,3)]), R, 1 ).^2, 2 ) );
    F_inc_l = @(R0,eps) Manifold.SO3.splus_l( R0, eps );
    F_inc_r = @(R0,eps) Manifold.SO3.splus_r( R0, eps );
    dist = deg2rad(90);
    gv = -dist:0.01:+dist;
    figure
    plotManifoldSubdim( gv, F_Cost, F_inc_l, Rgt );
    figure
    plotManifoldSubdim( gv, F_Cost, F_inc_l, Rtri );
    figure
    plotManifoldSubdim( gv, F_Cost, F_inc_r, Rgt );
    figure
    plotManifoldSubdim( gv, F_Cost, F_inc_r, Rtri );
end

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