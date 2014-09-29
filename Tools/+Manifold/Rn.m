classdef Rn < Manifold.Base
    %P2 Class for the Rn standard space (column vectors)
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = Rn( X, A_X )
            if nargin~=0
                m = size(X,2);
                N = size(X,1);
                if m == 1
                    obj.X = X;
                    obj.x = X;
                    
                    obj.dim = N;
                    obj.DIM = N;
                    
                    if exist('A_X','var')
                        obj.A_X = A_X;
                        obj.A_x = A_X;
                    end
                else
                    obj(1,m) = Manifold.Rn;
                    
                    X = num2cell( X, 1 );
                    [obj.X] = deal( X{:} );
                    [obj.x] = deal( X{:} );
                    
                    [obj.dim] = deal(N);
                    [obj.DIM] = deal(N);
                end
            end
        end
        
        function [out, ob] = plus( x1, x2 )
            if isa(x1,'Manifold.Rn')
                x1 = x1.X;
            end
            if isa(x2,'Manifold.Rn')
                x2 = x2.X;
            end
            [x1,x2] = Manifold.Rn.autodim( x1,x2 );
            out = x1 + x2;
            
            if nargout == 2
                ob = Manifold.Rn( out );
            end
        end
        function out = minus( x1, x2 )
            if isa(x1,'Manifold.Rn')
                x1 = x1.X;
            end
            if isa(x2,'Manifold.Rn')
                x2 = x2.X;
            end
            [x1,x2] = Manifold.Rn.autodim( x1,x2 );
            out = x1 - x2;
        end
        function out = mtimes( x1, x2 )
            if isa(x1,'Manifold.Rn')
                x1 = x1.X;
            end
            if isa(x2,'Manifold.Rn')
                x2 = x2.X;
            end
            out = x1 * x2;
        end
        
        function J_X_x = Dexp( obj )
            J_X_x = eye( obj.dim );
        end
        
        function J_x_X = Dlog( obj )
            J_x_X = eye( obj.dim );
        end
        
        function J_X_X = DLie( obj )
            J_X_X = eye( obj.dim );
        end
        
    end
    methods (Static)
        %         function n = exp( alpha )
        %             n = [ cos(alpha) ; sin(alpha) ];
        %         end
        %
        %         function alpha = log( n )
        %             alpha = atan2( n(2), n(1) );
        %         end
        %
        %         function mu_n = mean( X )
        %             mu_n = sum( X, 2 );
        %             mu_n = mu_n / norm( mu_n );
        %         end
        function mu = mean( in )
            if isa(in,'Manifold.Rn')
                in = [in.X];
            end
            mu = mean( in, 2 );
            mu = Manifold.Rn(mu);
        end
        
        function [x1,x2] = autodim( x1,x2 )
            % Auto repmat
            if size(x1,1)~=size(x2,1)
                error('Column vectors of different 1-dimension');
            end
            s1 = size(x1,2);
            s2 = size(x2,2);
            if s1 > 1 && s2 > 1 && s1~=s2
                error('Different non-single 2-dimensions');
            elseif s1 > 1 && s2 == 1
                x2 = repmat(x2,1,s1);
            elseif s1 == 1 && s2 > 1
                x1 = repmat(x1,1,s2);
            end
        end
    end
    
end