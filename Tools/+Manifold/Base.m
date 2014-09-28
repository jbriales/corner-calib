classdef Base < handle
    %CBASEMANIFOLD Abstract class for the features and methods in a
    %manifold
    %   Detailed explanation goes here
    
    properties
        X   % Manifold point
        x   % Minimal representation of point
        
        A_X % Covariance of value in representation space
        A_x % Minimal representation covariance (tangent space)
        
        dim % Dimension of the manifold (minimal representation)
        DIM % Dimension of the representation
    end
    
    methods % Not compulsory
        x = minus( obj, X )
        X = plus( obj, x )
        
        J = Dexp( obj )
        J = Dlog( obj )
        function J = Dproj( obj )
            % Useful result to constrain derivatives into manifold
            J = obj.Dexp * obj.Dlog;
        end
        
        function setMinimalCov( obj, A_x )
            obj.A_x = A_x;
            J = obj.Dexp;
            obj.A_X = J * A_x * J';
        end
        
        function setRepresentationCov( obj, A_X )
            obj.A_X = A_X;
            J = obj.Dlog;
            obj.A_x = J * A_X * J';
        end
    end
    
    methods (Static)
        % Static methods
        X = exp( x )
        x = log( X )
        m = mean( XX )
    end
end

