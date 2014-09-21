classdef S2 < Manifold.Base
    %S2 Class for the S2 (sphere) manifold
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = S2( n, A_n )
            if ~all( size(n) == [3 1] )
                error('Wrong inputs');
            end
            n = n / norm( n ); % Assure norm = 1
                        
            obj.X = n;
            obj.x = [];

            obj.dim = 2;
            obj.DIM = 3;
            
            if exist('A_n','var') && exist('A_p','var')
                obj.A_X = A_n;
                J_x_X = obj.Dlog;
                obj.A_x = J_x_X * A_n * J_x_X';
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