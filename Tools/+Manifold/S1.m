classdef S1 < Manifold.Base
    %S1 Class for the S1 (circle) manifold
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = S1( input )
            if all( size(input) == [2 1] )
                n = input / norm(input);
                a = atan2( n(2), n(1) );
            elseif numel( input ) == 1
                a = input;
                n = [ cos(a) ; sin(a) ];
            end
                
            obj.X = n;
            obj.x = a;
            obj.dim = 1;
            obj.DIM = 2;
        end
        
        function inc_alpha = minus( n, obj )
            n0 = obj.X;
            
            N = size( n, 2 );
            ort_n0 = [ -n0(2) +n0(1) ]';
            ort_n0 = repmat( ort_n0, 1, N );
            % TODO: asin or atan?
            inc_alpha = asin( dot( ort_n0, n, 1 ) );
        end
        
        function n = plus( obj, inc_alpha )
            N = size( inc_alpha, 2 );
            n  = zeros(2,N);
            n0 = obj.X;
            for i=1:N
                ca = cos(inc_alpha(i));
                sa = sin(inc_alpha(i));
                R = [ca -sa ; sa ca];
                n(:,i) = R * n0;
            end
        end
        
        function v = mtimes( M, obj )
            % Matrix multiplication of S1 vector
            v = M * obj.X;
        end
        
        function J_n_alpha = Dexp( obj )
            n = obj.X;
            J_n_alpha = [ -n(2) ; +n(1) ];
        end
        
        function J_alpha_n = Dlog( obj )
            n = obj.X;
            J_alpha_n = [ -n(2) , +n(1) ]; % Transpose of Dexp
        end
        
        function J_n_n = DLie( obj )
            % Useful result to constrain derivatives into manifold
            J_n_n = obj.Dexp * obj.Dlog;
        end
        
        % Static methods
    end
    methods (Static)        
        function n = exp( alpha )
            n = [ cos(alpha) ; sin(alpha) ];
        end
        
        function alpha = log( n )
            alpha = atan2( n(2), n(1) );
        end
        
        function mu_n = mean( X )
            mu_n = sum( X, 2 );
            mu_n = mu_n / norm( mu_n );
        end
    end
    
end
