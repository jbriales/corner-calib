classdef SO3 < Manifold.Base
    %CSO Class for the SO(3) manifold
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = SO3( R )
            obj.X = R;
            obj.x = obj.log( R );
            obj.dim = 3;
            obj.DIM = 9;
        end
        
        function inc_eps = minus( obj, X )
            R  = reshape( X, 3,3, [] );
            R0 = obj.X;
            
            N = size(X,2);
            inc_eps = zeros(3,N);
            for i=1:N
                inc_eps(:,i) = logmap( R(:,:,i) * R0' );
            end
        end
        
        function R = plus( obj, inc_eps )
            % inc_eps can be an array of column vectors
            N = size( inc_eps, 2 );
            R = zeros(3,3,N);
            R0 = obj.X;
            
            if size(inc_eps,2)>1 % Parallel computation
                if exist('multiprod','file') % manopt tool
                    R_eps = expmap( inc_eps );
                    R = multiprod( R_eps, R0 );
                else
                    for i=1:N
                        R(:,:,i) = expmap( inc_eps(:,i) ) * R0;
                    end
                end
            else
                R = expmap( inc_eps ) * R0;
            end
        end
        
        function J_R_eps = Dexp( obj )
            R = obj.X;
            J_R_eps  = [ -skew(R(:,1))
                         -skew(R(:,2))
                         -skew(R(:,3)) ];
        end
        
        function J_eps_R = Dlog( obj )
            R = obj.X;
            J_eps_R = 0.5 * [ -skew(R(:,1))
                              -skew(R(:,2))
                              -skew(R(:,3)) ]';
        end
        
        % Static methods
    end
    methods (Static)        
        function R = exp( u )
            %  R = exp( inc_eps )
            %  Rodrigues' rotation formula.
            if nargin==1
                theta = norm(u,2);
            end
            
            if theta == 0
                R = eye(3);
            else
                u = u./norm(u,2);
                S = [    0 -u(3)  u(2);
                    u(3)   0  -u(1);
                    -u(2) u(1)   0  ];
                R = eye(3) + sin(theta)*S + (1-cos(theta))*S^2;
            end
        end
        
        function u = log( R )
            %  inc_eps = log( R )
            %  The Rodrigues' formula for rotation matrices.
            [~,~,V] = svd(R - eye(3));
            u = V(:,3);
            
            c2 = trace(R) - 1; % cos(theta)*2
            s2 = u' * [R(3,2)-R(2,3) , R(1,3)-R(3,1) , R(2,1)-R(1,2)]';
            theta = atan2(s2,c2);
            
            u = u * theta;
        end
        
        function mu_R = mean( RR )
            R  = reshape( RR, 3,3, [] );
            [U,~,V] = svd( sum( R, 3 ) );
            mu_R = U*V';
        end
        
        function R = splus_l( R0, inc_eps )
            % inc_eps can be an array of column vectors
            N = size( inc_eps, 2 );
            R = zeros(3,3,N);
            
            if size(inc_eps,2)>1 % Parallel computation
                if exist('multiprod','file') % manopt tool
                    R_eps = expmap( inc_eps );
                    R = multiprod( R_eps, R0 );
                else
                    for i=1:N
                        R(:,:,i) = expmap( inc_eps(:,i) ) * R0;
                    end
                end
            else
                R = expmap( inc_eps ) * R0;
            end
        end
        function R = splus_r( R0, inc_eps )
            % inc_eps can be an array of column vectors
            N = size( inc_eps, 2 );
            R = zeros(3,3,N);
            
            if size(inc_eps,2)>1 % Parallel computation
                if exist('multiprod','file') % manopt tool
                    R_eps = expmap( inc_eps );
                    R = multiprod( R0, R_eps );
                else
                    for i=1:N
                        R(:,:,i) = expmap( inc_eps(:,i) ) * R0;
                    end
                end
            else
                R = expmap( inc_eps ) * R0;
            end
        end
    end
    
end

