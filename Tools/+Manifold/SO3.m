classdef SO3 < Manifold.Base
    %CSO Class for the SO(3) manifold
    %   Detailed explanation goes here
    
    properties
        R   % Property to store 3x3 matrix
    end
    
    methods
        function obj = SO3( R )
            if nargin ~= 0 % To validate no input call
                m = size(R,3); % Number of inputs in 3rd dimension
                if m == 1
                    % This case is maintaned for the sake of clarity
                    % TODO: X should be only column vector for uniformity with other manifolds
                    obj.X = R;
                    obj.x = obj.log( R );
                    obj.dim = 3;
                    obj.DIM = 9;
                else
                    % Parallel constructor case
                    obj(1,m) = Manifold.SO3;
                    
                    % Compute auxiliary values
                    u = obj.log( R );
                    
                    % Create cell of 3x3 Rs, removing singleton dimensions
                    % This type allows to work with lists
                    R = squeeze( num2cell( R, [1 2] ) );
                    u = num2cell( u, 1 );

                    % Assign
                    [obj.X] = deal( R{:} );
                    [obj.x] = deal( u{:} );
                    
                    [obj.dim] = deal(3);
                    [obj.DIM] = deal(9);
                end
            end
        end
        
        function eps = minus( R1, R2 )
            if isa(R1,'Manifold.SO3')
                R1 = R1.X;
            end
            if isa(R2,'Manifold.SO3')
                R2 = R2.X;
            end
            eps = Manifold.SO3.sminus( R1, R2 );
        end
        
        function [R, obj_out] = plus( obj, inc_eps )
            R0 = obj.X;
            R = obj.splus_l( R0, inc_eps );
            
            if nargout == 2
                obj_out = Manifold.SO3( R );
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
        
        function J_R_R = Dproj( obj )
            J_R_R = obj.Dexp * obj.Dlog;
        end
        
        % Static methods
    end
    methods (Static)        
        function R = exp( u )
            % R = exp( u )
            % The Rodrigues' formula for rotation matrices.
            % u can be an array of column vectors
            R = expmap( u );
        end
        
        function u = log( R )
            % u = log( R )
            % The Rodrigues' formula for rotation matrices.
            % R can be an 3x3xN array of rotation matrices
            u = logmap( R );
        end
        
        % NOTE: Static methods should substitute to class method as factory
        % of methods
        function J_R_eps = Jexp( R )
            % J_R_eps = Jexp( R )
            % Returns 9x3 jacobian of Matrix SO(3) wrt minimal so(3)
            % representation
            J_R_eps  = [ -skew(R(:,1))
                         -skew(R(:,2))
                         -skew(R(:,3)) ];
        end
        
        % NOTE: Static methods should substitute to class method as factory
        % of methods
        function J_eps_R = Jlog( R )
            % J_R_eps = Jlog( R )
            % Returns 3x9 jacobian of minimal so(3) representation wrt Matrix SO(3)
            J_eps_R = 0.5 * [ -skew(R(:,1))
                              -skew(R(:,2))
                              -skew(R(:,3)) ]';
        end
        
        function mu_R = mean( RR )
            if isa(RR,'Manifold.SO3')
                RR = [RR.X];
            end
            R  = reshape( RR, 3,3, [] );
            [U,~,V] = svd( sum( R, 3 ) );
            mu_R = U*V';
            mu_R = Manifold.SO3( mu_R );
        end
        
        function eps = sminus( R, R0 )
            if exist('multiprod','file') % manopt tool
                R_eps = multiprod( R, multitransp(R0) );
                eps   = Manifold.SO3.log( R_eps );
            else
                N = size(R,2);
                eps = zeros(3,N);
                for i=1:N
                    eps(:,i) = logmap( R(:,:,i) * R0' );
                end
            end
        end
        
        function R = splus_l( R0, inc_eps )
            % inc_eps can be an array of column vectors
            N = size( inc_eps, 2 );
            R = zeros(3,3,N);
            
            R_eps = expmap( inc_eps );
            if size(inc_eps,2)>1 % Parallel computation
                if exist('multiprod','file') % manopt tool
                    R = multiprod( R_eps, R0 );
                else
                    for i=1:N
                        R(:,:,i) = R_eps(:,:,i) * R0;
                    end
                end
            else
                R = R_eps * R0;
            end
        end
        function R = splus_r( R0, inc_eps )
            % inc_eps can be an array of column vectors
            N = size( inc_eps, 2 );
            R = zeros(3,3,N);
            
            R_eps = expmap( inc_eps );
            if size(inc_eps,2)>1 % Parallel computation
                if exist('multiprod','file') % manopt tool
                    R = multiprod( R0, R_eps );
                else
                    for i=1:N
                        R(:,:,i) = R0 * R_eps(:,:,i);
                    end
                end
            else
                R = R0 * R_eps;
            end
        end
    end
    
end

