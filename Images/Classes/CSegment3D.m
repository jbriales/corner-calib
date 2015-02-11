classdef CSegment3D
    %CSegment3D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        p1
        p2
        v
    end
    
    methods
        function this = CSegment3D(p1,p2)
            this.p1 = p1(:);
            this.p2 = p2(:);
            
            this.v  = snormalize(p2-p1);
        end
        
        function seg = inverse( this )
            seg = CSegment3D(this.p2, this.p1);
        end
        
        function transform( this, type, T )
            switch type
                case 'euclidean3D'
                    R = T(:,1:3);
                    t = T(:,4);
                    this.p1 = R * this.p1  + t;
                    this.p2 = R * this.p2  + t;
                otherwise
                    warning('Unknown transformation %s\n',type);
            end
        end
        
        function seg2D = project( this, Cam )
            % Array-aware operation 
            if numel(this)==1
                pixels = Cam.projectPoints( [this.p1, this.p2] );
                seg2D = CSegment2D( pixels(:,1),pixels(:,2) );
            else
                pixels = Cam.projectPoints( [this.p1, this.p2] );
                C = mat2cell(pixels, 2, [numel(this) numel(this)]);
                seg2D = CSegment2D( C{:} );
            end
        end
        
        function v = projectVector( this, Cam )
            pixels = Cam.projectPoints( [this.p1, this.p2] );
            v = snormalize( pixels(:,2) - pixels(:,1) );
        end
        
        function l = projectLine( this, Cam )
            if numel(this)>1
                n = numel(this);
                cell_l = cell(1,n);
                for k=1:n
                    cell_l{k} = this(k).projectLine(Cam);
                end
                l = cell_l; % Return cell array or array?
                return % Finish plot
            end
                
            pixels = Cam.projectPoints( [this.p1, this.p2] );
            l = cross( makehomogeneous(pixels(:,1)),...
                       makehomogeneous(pixels(:,2)) );
        end
        
        % Overload method for plot=plot3
        function [h,g] = plot( this, format, tag )
            [h,g] = this.plot3( format, tag );
        end
        
        function [h,g] = plot3( this, format, tag )
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
                        h(k) = this(k).plot3(format);
                    else
                        [h(k),g(k)] = this(k).plot3(format,tag{k});
                    end
                end
                return % Finish plot
            end
            
            vec = [this.p1(:) this.p2(:)]';
            cell_vec = num2cell(vec,[1,3]);
            if iscell( format )
                h = plot3( cell_vec{:}, format{:} );
            else
                h = plot3( cell_vec{:}, format );
            end
            if ~isempty(tag)
                xyz = 0.5*(this.p1 + this.p2);
                g = text(xyz(1),xyz(2),xyz(3),tag);
            end
        end
    end
    
end

