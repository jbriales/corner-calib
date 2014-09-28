classdef Cxi < Manifold.Dyn
    %Cxi Class for storage of Trihedron Observation image, structure formed
    %by one point (center) and three directions (one per line)
    %   Detailed explanation goes here
    
    properties (Dependent)
%         c   % 2x1 R2 vector
%         v   % 1x3 array with S1 directions of lines
    end
    
    methods % Not compulsory
        function obj = Cxi( c, v1, v2, v3 )
            obj = obj@Manifold.Dyn( c, v1, v2, v3 );
        end
        
        % Function to simplify interface
        function c = c( obj )
            c = obj.vars{1}.X;
        end
        function v = v( obj, ind )
            v = obj.vars{1+ind};
        end
    end
end
