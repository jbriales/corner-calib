classdef CSegment2D < handle
    %CSegment3D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        p1
        p2
        
        % Covariance
        A_l % Covariance of homogeneous lines (3x3 rank-defficient matrix)
        
        % Auxiliar
        tag % Allow to tag a segment object in order to identify it in matching processes, etc.
    end
    
    properties(Dependent)
        c % Central point (mean of p1 and p2)
        v % Direction vector
        n % Normal vector
        l % Homogeneous line
        d % Segment length
    end
    
    methods
        function this = CSegment2D(p1,p2)           
            if nargin ~= 0 % Allow nargin == 0 syntax
                if any( size(p1) ~= size(p2) ) || size(p1,1) ~= 2
                    error('Bad input dimensions');
                end
                
                n = size(p1,2);
                this(n) = CSegment2D; % Preallocate object array
                for j = 1:n
                    this(j).p1 = p1(:,j);
                    this(j).p2 = p2(:,j);
                end
            end
        end
        
        % Overloaded operations
        function p = mtimes( this, other )
            assert( strcmp( class(this), class(other) ) );
            if numel(other) == 1
                p = makeinhomogeneous( cross(this.l,other.l) );
            else
                p = makeinhomogeneous( skew(this.l) * [other.l] );
            end
        end
        function d = minus( pts, this )
            assert( size(pts,1)==2 );
            m_ = [this.v; -this.v'*this.c];
            d = m_' * makehomogeneous( pts );
        end
        
        function transform( this, type, T )
            switch type
                case 'euclidean2D'
                    R = T(:,1:2);
                    t = T(:,3);
                    this.p1 = R * this.p1  + t;
                    this.p2 = R * this.p2  + t;
                otherwise
                    warning('Unknown transformation %s\n',type);
            end
        end
        
        % Different segment properties (obtained from p1 and p2)
        function c = get.c(this)
            c = 0.5*(this.p1+this.p2);
        end
        function v = get.v(this)
            v = snormalize( this.p2-this.p1 );
        end
        function n = get.n(this)
            n = [ 0 -1 ; 1 0 ] * this.v;
        end
        function d = get.d(this)
            d = norm( this.p2 - this.p1 );
        end
        function l = get.l(this)
            l = [ this.n; -this.n'*this.p1 ];
        end
        
        function A_l = get.A_l( this )
            % A_l = get.A_l( this )
            % The A_l property is only computed once and then stored to
            % A_l property of object.
            % If A_l is called after computation, it...
            
            if isempty(this.A_l)
                % Uncertainty computation based on segment length inverse
                Jac_l_v = [ eye(2)
                            -this.c' ] * this.v; % Projected to circle
                A_v = 1 / this.d; % Ficticial, make uncertainty inv prop to length
                A_l = Jac_l_v * A_v * Jac_l_v';
                % NOTE: Center covariance has been disregarded here
                this.A_l = A_l;
            else
                A_l = this.A_l;
            end
        end
            
        function A_N = estimateNcovLength( this, segs )
            % A_N = estimateNcovLength( this, segs )
            % Get uncertainty estimation from segmens and length
            
            % Produce covariance of segs wrt its length
            A_v = cell(1,3);
            A_l = cell(1,3);
            A_n = cell(1,3);
            for k=1:3
                Jac_l_v = [ eye(2)
                            -segs(k).c' ] * segs(k).v; % Projected to circle
                A_v{k} = 1 / segs(k).d; % Ficticial, make uncertainty inv prop to length
                A_l{k} = Jac_l_v * A_v{k} * Jac_l_v';
                % NOTE: Center covariance has been disregarded here
                
                % Compute normal covariance from line covariance
                Jac_n_l = Dsnormalize(this.K'*segs(k).l) * this.K';
                A_n{k} = Jac_n_l * A_l{k} * Jac_n_l';
            end
            A_N = blkdiag( A_n{:} );
        end

        function [h,g] = plot( this, format, tag )
            % [h,g] = this.plot( format, tag )
            % Receives:
            %       format - for segment lines (as in plot)
            %       tag - string (or cell of strings for array object)
            %   with values to print next to segments
            % Returns:
            %       h - array of handles for printed lines
            %       g - array of handles for added tags            
            
            % Complete input with default values
            if ~exist('format','var')
                format = '-k';
            end
            if ~exist('tag','var')
                tag = [];
            end
            
            if numel(this)>1
                n = numel(this);
                % Solve multiple cases with recursive calls
                h = zeros(1,n);
                g = zeros(1,n);
                for k=1:n
                    if isempty(tag)
                        h(k) = this(k).plot(format);
                    else
                        if islogical( tag )
                            % Option true (logical)
                            [h(k),g(k)] = this(k).plot(format,true);
                        elseif iscell( tag )
                            [h(k),g(k)] = this(k).plot(format,tag{k});
                        else
                            warning('Check tag option');
                            keyboard;
                        end
                    end
                end
                return % Finish plot
            end
            
            vec = [this.p1(:) this.p2(:)]';
            cell_vec = num2cell(vec,[1,3]);
            if ~iscell(format)
                % To use common interface for format options
                format = {format};
            end
            h = plot( cell_vec{:}, format{:} );
            
            %% Plot segment tag
            % Pretreat tag data (use given input or object property)
            if tag == true
                tag = this.tag;
                if isnumeric(tag)
                    tag = num2str(tag);
                end
            end
            g = NaN; % Default value if tag is not plot
            if ~isempty(tag)
                space = 0.01;
                xy = 0.5*(this.p1 + this.p2) + space*this.n;
                g = text(xy(1),xy(2),tag);
            end
        end
    end
    
end