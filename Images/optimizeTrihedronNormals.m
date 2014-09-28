function obj_Rtri = optimizeTrihedronNormals( obj_Nbp, labels, R0 )
% obj_Rtri = otpimizeTrihedronNormals( obj_Nbp, R0 )
% Function that receives a Dyn Manifold with a set of back-projected plane
% normals and an initial estimate (which can be obtained from the minimal
% solution computeTrihedronNormals)
% An array of labels (1,2,3) is also given to tag the world direction of
% imaged lines

N = obj_Nbp.arr;
canonical = zeros(3,size(N,2)); % Array with canonical direction corresponding to each line
canonical(1,labels==1) = 1;
canonical(2,labels==2) = 1;
canonical(3,labels==3) = 1;
F_Phi = @(R) dot( N, R*canonical, 1 )';
J_Phi = @(R) cross(  R*canonical, N, 1 )';

[ Rtri, err, errNorm, W ] = LM_Man_optim(@LM_Fun,R0,...
            'space','SO(3)','weighted',true);%,...
%             'debug',obj.debug_level,...
%             'maxIters', obj.maxIters,...
%             'minChange', obj.minParamChange,...
%             'minErrorChange', obj.minErrorChange);
obj_Rtri = Manifold.SO3( Rtri );

% Compute covariance of trihedron normals
% The covariance keeps the same even if more obs are added? Error?
J_Phi_eps = J_Phi( Rtri );

J_Phi_Nbp = num2cell( Rtri*canonical, 1 );
J_Phi_Nbp = blkdiag( J_Phi_Nbp{:} )';

J_eps_Nbp = - J_Phi_eps \ J_Phi_Nbp;
% J_Rtri_eps = obj_Rtri.Dexp;
% J_Rtri_Nbp = J_Rtri_eps * ( - J_phi_eps \ J_Phi_Nbp );

obj_Rtri.setMinimalCov( J_eps_Nbp * obj_Nbp.A_X * J_eps_Nbp' );

    function [res, jac, weights] = LM_Fun( R )
%         N_to = size( Nbp, 2 ) / 3;
%         N_tri = R * kron(ones(1,N_to),eye(3));
%         
%         res = dot( Nbp, N_tri, 1 )';
%         jac = cross( N_tri, Nbp )';
        res = F_Phi( R );
        jac = J_Phi( R );
        weights = eye(length(res));
    end
end

