classdef Dyn < Manifold.Base
    %CDyn Class for storage of cartesian product manifold, structure formed
    %by several manifold variables
    %   Detailed explanation goes here
    
    properties
        vars    % Variables composing the manifold
        Nvars   % Number of variables in manifold product
        
        idxs    % Cell array with indexes corresponding to minimal repr.
        IDXS    % Cell array with indexes corresponding to complete repr.
    end
    
    methods % Not compulsory
        function obj = Dyn( varargin )
            if nargin ~= 0
                obj.X = zeros(0,1);
                obj.dim = 0;
                obj.DIM = 0;
                count = 1;
                for i=1:nargin
                    % If any input is an object array, take element wise
                    in_objs = varargin{i};
                    for k=1:numel(in_objs)
                        ob = in_objs(k);
                        obj.vars{count} = ob;
                        obj.X = [ obj.X ; ob.X(:) ];
                        
                        obj.idxs{count} = obj.dim + (1:ob.dim);
                        obj.IDXS{count} = obj.DIM + (1:ob.DIM);
                        
                        obj.dim = obj.dim + ob.dim;
                        obj.DIM = obj.DIM + ob.DIM;
                        
                        count = count + 1; % Increase total counter
                    end
                end
                obj.Nvars = numel( obj.vars );
            end
        end
        
        function obj = setFromX( obj, X )
            % obj = setFromX( obj, X )
            % Method to set Dyn manifold using existing reference with
            % input a set of vectorized values
            % TODO
            for k=1:obj.Nvars
                
            end
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
        
        function [out_X,out_ob] = plus( obj, inc_eps )
            N_eps = size(inc_eps,2);
            temp = cell(1,obj.Nvars);
            out_X = [];
            for k=1:obj.Nvars
                X = obj.vars{k} + inc_eps(obj.idxs{k},:);
                out_X = [out_X ; X];
                eval(strcat('ob=',class(obj.vars{k}),'(X);'));
                temp{k} = ob;
            end

            if nargout == 2
%                 out_ob(1,N_eps) = Manifold.Dyn;
%                 eval( strcat('out_ob(1,N_eps) = ',class(obj),';') );
                for i=1:N_eps
                    c = cell(1,obj.Nvars);
                    for j=1:obj.Nvars
                        c{j} = temp{j}(i);
                    end
                    eval( strcat('out_ob(i) = ',class(obj),'(c{:});') );
%                     out_ob(i) = Manifold.Dyn(c{:});
                    % Set same cov as input to avoid errors in input functions
                    if ~isempty( obj.A_x ) % This should be done elsewhere
                        out_ob(i).setMinimalCov(obj.A_x);
                    end
                end
            end
        end
        
        function inc_eps = minus( obj, in )
            all_inc = [];
            for k=1:obj.Nvars
                inc = obj.vars{k} - in.vars{k};
                all_inc = [all_inc ; inc];
            end
            inc_eps = all_inc;
        end
        
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
    
    methods (Static)
        % Static methods
%         X = exp( x )
%         x = log( X )
        function out = mean( objs )
            vars = cell(1,objs(1).Nvars);
            % Extract different manifolds
            temp = reshape( [ objs.vars ], objs(1).Nvars,[] );
            for k=1:objs(1).Nvars
                eval(strcat('var_mu=',...
                     class(objs(1).vars{k}),'.mean([temp{k,:}])'));
                vars{k} = var_mu;
            end
            eval(strcat('out=',class(objs),'(vars{:});'));
%             out = Manifold.Dyn( vars{:} );
        end
    end
end