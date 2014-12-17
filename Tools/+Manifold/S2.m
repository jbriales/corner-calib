classdef S2 < Manifold.Base
    %S2 Class for the S2 (sphere) manifold
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = S2( n, A_n )
            if nargin ~= 0
                m = size(n,2);
                if m == 1
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
                else
                    obj(1,m) = Manifold.S2;
                    
                    n = num2cell( n, 1 );
                    
                    % Assign
                    [obj.X] = deal( n{:} );
                    
                    [obj.dim] = deal(2);
                    [obj.DIM] = deal(3);
                end
            end
        end
        
        function [n, ob] = plus( obj, eps )
            H = Householder( obj.X );
            % Check H projects to [0,0,+1]
            % Is correct to multiply by -1?
            if [0 0 1] * H * obj.X < 0
                H = - H;
            end
            n_eps = obj.exp( eps );
            
            n = H * n_eps;            
            
            if nargout == 2
                ob = Manifold.S2( n );
            end
        end
        
        function eps = minus( n1, n2 )
            if isa(n1,'Manifold.S2')
                n1 = n1.X;
            end
            if isa(n2,'Manifold.S2')
                n2 = n2.X;
            end
            n  = n1; % The element being substracted
            n0 = n2; % The element substracting
            
            H = Householder( n0 );
            n_eps = H * n;
            eps = Manifold.S2.log( n_eps );
        end
        
        % Updated to Householder
        function J_X_x = Dexp( obj )
            H = Householder( obj.X );
            J_X_x = H * [eye(2) zeros(2,1)]';
        end
        
        function J_x_X = Dlog( obj )
            H = Householder( obj.X );
            J_x_X = [eye(2) zeros(2,1)] * H;
        end
        
        function J_X_X = DLie( obj )
            J_X_X = obj.Dexp * obj.Dlog;
        end
%         function J_X_x = Dexp( obj )
%             l = obj.X;
%             J_X_x = null( l' ); % Transpose of Dlog
%         end
%         
%         function J_x_X = Dlog( obj )
%             l = obj.X;
%             J_x_X = null( l' )';
%         end
%         
%         function J_X_X = DLie( obj )
%             J_X_X = obj.Dexp * obj.Dlog;
%         end
        
    end
    methods (Static)        
        function n = exp( eps )
            % n = exp( eps )
            % Uses exp map explained in Zisserman
            eps_ = makehomogeneous( eps );
            n = snormalize(eps_);
        end
        
        function eps = log( n )
            eps = makeinhomogeneous( n );
        end
        
        function mu_n = mean( in )
            if isa(in,'Manifold.S2')
                in = [in.X];
            end
            mu_n = sum( in, 2 );
            mu_n = mu_n / norm( mu_n );
            mu_n = Manifold.S2( mu_n );
        end
    end
    
end