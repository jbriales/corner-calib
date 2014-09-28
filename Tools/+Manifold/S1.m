classdef S1 < Manifold.Base
    %S1 Class for the S1 (circle) manifold
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = S1( input )
            if nargin ~= 0 % To validate no input call
                m = size(input,2); % Number of inputs in 3rd dimension
                if m == 1
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
                else
                    if size(input,1) == 2
                        n = snormalize( input );
                        a = atan2( n(2,:,:), n(1,:,:) );
                    else % Make more robust the access to this option
                        a = input;
                        n = [ cos(a) ; sin(a) ];
                    end
                    % Parallel constructor case
                    obj(1,m) = Manifold.S1;
                    
                    n = num2cell( n, 1 );
                    a = num2cell( a );
                    % Assign
                    [obj.X] = deal( n{:} );
                    [obj.x] = deal( a{:} );
                    [obj.dim] = deal(1);
                    [obj.DIM] = deal(2);
                end
            end
        end
        
        function inc_alpha = minus( n1, n2 )
            if isa(n1,'Manifold.S1')
                n1 = n1.X;
            end
            if isa(n2,'Manifold.S1')
                n2 = n2.X;
            end
            n  = n1; % The element being substracted
            n0 = n2; % The element substracting
            
            N = size( n1, 2 );
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
        
        % Static methods
    end
    methods (Static)        
        function n = exp( alpha )
            n = [ cos(alpha) ; sin(alpha) ];
        end
        
        function alpha = log( n )
            alpha = atan2( n(2), n(1) );
        end
        
        function mu_n = mean( in )
            if isa(in,'Manifold.S1')
                in = [in.X];
            end
            mu_n = sum( in, 2 );
            mu_n = mu_n / norm( mu_n );
            mu_n = Manifold.S1( mu_n );
        end
    end
    
end
