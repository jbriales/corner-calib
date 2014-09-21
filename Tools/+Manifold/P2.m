classdef P2 < Manifold.Base
    %P2 Class for the P2 (projective) manifold
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = P2( in1, p, A_n, A_p )
            if nargin == 1 % Input l
                l = in1;
            else
                n = in1;
                if ~all( size(n) == [2 1] & size(p) == [2 1] )
                    error('Wrong inputs');
                end
                l = [ n ; -n*p ];
            end
                        
            obj.X = l;
            obj.x = [];

            obj.dim = 2;
            obj.DIM = 3;
            
            if exist('A_n','var') && exist('A_p','var')
                J = [ eye(2) , zeros(2,2);
                    -p'    , -n'       ];
                A_l = J * blkdiag( A_n, A_p ) * J';
                
                obj.A_X = A_l;
                J_x_X = obj.Dlog;
                obj.A_x = J_x_X * A_l * J_x_X';
            end
        end
        
        
        function J_X_x = Dexp( obj )
            l = obj.X;
            J_X_x = null( l' ); % Transpose of Dlog
        end
        
        function J_x_X = Dlog( obj )
            l = obj.X;
            J_x_X = null( l' )';
        end
        
        function J_X_X = DLie( obj )
            J_X_X = obj.Dexp * obj.Dlog;
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
    end
    
end