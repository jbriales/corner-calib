classdef Csegment < Manifold.Dyn
    %Cxi Class for storage of Trihedron Observation image, structure formed
    %by one point (center) and three directions (one per line)
    %   Detailed explanation goes here
    
    properties (Dependent)
%         c   % 2x1 R2 vector
%         v   % 1x3 array with S1 directions of lines
    end
    
    methods % Not compulsory
        function obj = Csegment( c, v )
            % obj = Cxi( c, v1, v2, v3 )
            if isnumeric(c)
                c = Manifold.Rn( c );
            end
            if isnumeric(v)
                v = Manifold.Rn( v );
            end
            
            % obj = Cxi( c_R2, v1_S1, v2_S1, v3_S1 )
            obj = obj@Manifold.Dyn( c, v );
        end
        
        % Function to simplify interface
        function c = c( obj )
            c = obj.vars{1}.X;
        end
        function vec = vec( obj )
            vec = obj.vars{2}.X;
        end
        function v = v( obj )
            vec = obj.vars{2}.X;
            v = snormalize(vec);
        end
        function r = r( obj )
            vec = obj.vars{2}.X;
            r = norm(vec);
        end
    end
end
