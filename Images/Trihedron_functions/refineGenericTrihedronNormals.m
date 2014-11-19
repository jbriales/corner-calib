function rho = refineGenericTrihedronNormals( D, c, rho0, debug )
% rho = refineGenericTrihedronNormals( D, c, rho0 )
% Numerically solve the system of equations for trihedron directions
% when generic angles exist between lines (assumed near orthogonality)

if ~exist('debug','var')
    debug = true;
end

rho0 = mat2cell( rho0, 2,[1 1 1] );

null2 = @(x)[-x(2); x(1)]; % Compute orthogonal vector in positive orientation
function [F,JF] = Fun( rho )
    rho = mat2cell( rho, 2,[1 1 1] );
    F  = [ rho{1}'*D{3}*rho{2} - c(3)
           rho{2}'*D{1}*rho{3} - c(1)
           rho{3}'*D{2}*rho{1} - c(2) ];
    JF = [ rho{2}'*D{3}*null2(rho{1}') , rho{1}'*D{3}*null2(rho{2}') , 0
           0 , rho{3}'*D{1}*null2(rho{2}') , rho{2}'*D{1}*null2(rho{3}')
           rho{3}'*D{2}*null2(rho{1}') , 0 , rho{1}'*D{2}*null2(rho{3}') ];
end

obj_rho = Manifold.Dyn( Manifold.S1(rho0{1}),...
                        Manifold.S1(rho0{2}),...
                        Manifold.S1(rho0{3}) );
                    
[ obj_rho, Err ] = Newton_solve( @(var)Fun(reshape(var.X,2,3)),...
    obj_rho, 'space', 'Manifold', 'debug', debug );
 rho = reshape( obj_rho.X, 2,3 );
end