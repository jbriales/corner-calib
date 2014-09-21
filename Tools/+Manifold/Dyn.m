classdef Dyn < Manifold.Base
    %CDyn Class for storage of dynamic manifold, structure formed
    %by several manifold variables
    %   Detailed explanation goes here
    
    properties
        c   % 2x1 R2 vector
        v   % 1x3 array with S1 directions of lines
        
        vars    % Variables composing the manifold
        Nvars   % Number of variables in composed manifold
    end
    
    methods % Not compulsory
        function obj = Dyn( varargin )
            obj.X = zeros(0,1);
            obj.dim = 0;
            obj.DIM = 0;
            for i=1:nargin
                obj.vars{i} = varargin{i};
                obj.X = [ obj.X ; varargin{i}.X(:) ];
                obj.dim = obj.dim + varargin{i}.dim;
                obj.DIM = obj.DIM + varargin{i}.DIM;
            end
            obj.Nvars = nargin;
        end
        
        function arr = arr( obj )
            % Return composed manifold variables as an array of columns (if
            % possible)
            N = obj.Nvars;
            DIMS = zeros(1,N);
            arr  = cell(1,N);
            for i=1:N
                DIMS(i) = obj.vars{i}.DIM;
                arr{i}  = obj.vars{i}.X;
            end
            if all( repmat(DIMS(1),1,N) == DIMS )
                arr = cell2mat( arr );
            end
        end
        
%         x = minus( obj, X )
%         X = plus( obj, x )
        
        function J = Dexp( obj )
            J = cell(1,obj.Nvars);
            for i=1:obj.Nvars
                J{i} = obj.vars{i}.Dexp;
            end
            J = blkdiag( J{:} );
        end        
        function J = Dlog( obj )
            J = cell(1,obj.Nvars);
            for i=1:obj.Nvars
                J{i} = obj.vars{i}.Dlog;
            end
            J = blkdiag( J{:} );
        end
        function J = DLie( obj )
            J = obj.Dexp * obj.Dlog;
        end
    end
    
%     methods (Static)
%         % Static methods
%         X = exp( x )
%         x = log( X )
%         m = mean( XX )
%     end
end