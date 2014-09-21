classdef Cxi < Manifold.Base
    %Cxi Class for storage of Trihedron Observation image, structure formed
    %by one point (center) and three directions (one per line)
    %   Detailed explanation goes here
    
    properties
        c   % 2x1 R2 vector
        v   % 1x3 array with S1 directions of lines
    end
    
    methods % Not compulsory
        function obj = Cxi( c, v1, v2, v3 )
            obj.c = c;
            
            obj.v = Manifold.S1.empty(1,0);
            dirs = {v1, v2, v3};
            for i = 1:3
                if isa(dirs{i},'Manifold.S1')
                    obj.v(i) = dirs{i};
                else
                    obj.v(i) = Manifold.S1( dirs{i} );
                end
            end
            
            obj.X = [ obj.c , obj.v(1).X, obj.v(2).X, obj.v(3).X ];
            obj.X = obj.X(:);
            
            %
            
            obj.dim = 5;
            obj.DIM = 8;
        end      
        
%         x = minus( obj, X )
%         X = plus( obj, x )
        
        function J = Dexp( obj )
            J = blkdiag( eye(2), obj.v(1).Dexp, obj.v(2).Dexp, obj.v(3).Dexp );
        end
        function J = Dlog( obj )
            J = blkdiag( eye(2), obj.v(1).Dlog, obj.v(2).Dlog, obj.v(3).Dlog );
        end
%         J = Dlog( obj )
    end
    
    methods (Static)
        % Static methods
        X = exp( x )
        x = log( X )
        m = mean( XX )
    end
end
