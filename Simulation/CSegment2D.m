classdef CSegment2D
    %CSegment3D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        p1
        p2
        v
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
        
        function h = plot( this, format )
            % Complete input with default values
            if ~exist('format','var')
                format = '-k';
            end
            
            if numel(this)>1
                n = numel(this);
                % Solve multiple cases with recursive calls
                h = zeros(1,n);
                for k=1:n
                    h(k) = this(k).plot(format);
                end
                return % Finish plot
            end
            
            vec = [this.p1(:) this.p2(:)]';
            cell_vec = num2cell(vec,[1,3]);
            h = plot( cell_vec{:}, format );
        end
    end
    
end