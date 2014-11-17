function rho = refineGenericTrihedronNormals2( D, c, rho0 )
% rho = refineGenericTrihedronNormals( D, c, rho0 )
% Numerically solve the system of equations for trihedron directions
% when generic angles exist between lines (assumed near orthogonality)

rho0 = mat2cell( rho0, 2,[1 1 1] );

function [F2,JF2] = Fun( rho )
    rho = mat2cell( rho, 2,[1 1 1] );
    F  = [ rho{1}'*D{3}*rho{2} - c(3)
           rho{2}'*D{1}*rho{3} - c(1)
           rho{3}'*D{2}*rho{1} - c(2) ];
    F2 = 0.5*F'*F;
    JF = [ rho{2}'*D{3}*null(rho{1}') , rho{1}'*D{3}*null(rho{2}') , 0
           0 , rho{3}'*D{1}*null(rho{2}') , rho{2}'*D{1}*null(rho{3}')
           rho{3}'*D{2}*null(rho{1}') , 0 , rho{1}'*D{2}*null(rho{3}') ];
    JF2 = F'*JF;
end

obj_rho = Manifold.Dyn( Manifold.S1(rho0{1}),...
                        Manifold.S1(rho0{2}),...
                        Manifold.S1(rho0{3}) );
                    
[ obj_rho, Err ] = LM_Man_optim( @(var)Fun(reshape(var.X,2,3)),...
    obj_rho, 'space', 'Manifold' );
 rho = reshape( obj_rho.X, 2,3 );
end