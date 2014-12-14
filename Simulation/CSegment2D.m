classdef CSegment2D
    %CSegment3D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        p1
        p2
        v
    end
    
    properties(Dependent)
        n
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
                    this(j).v  = snormalize( p2(:,j)-p1(:,j) );
                end
            end
        end
        
        function n = get.n(this)
            n = [ 0 -1 ; 1 0 ] * this.v;
        end
        
        function [h,g] = plot( this, format, tag )
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
                        [h(k),g(k)] = this(k).plot(format,tag{k});
                    end
                end
                return % Finish plot
            end
            
            vec = [this.p1(:) this.p2(:)]';
            cell_vec = num2cell(vec,[1,3]);
            h = plot( cell_vec{:}, format );
            if ~isempty(tag)
                xy = 0.5*(this.p1 + this.p2);
                g = text(xy(1),xy(2),tag);
            end
        end
    end
    
end