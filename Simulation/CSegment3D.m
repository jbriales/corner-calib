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
            pixels = Cam.projectPoints( [this.p1, this.p2] );
            l = cross( makehomogeneous(pixels(:,1)),...
                       makehomogeneous(pixels(:,2)) );
        end
        
        function h = plot3( this, format )
            if ~exist('format','var')
                format = '-k';
            end
            vec = [this.p1(:) this.p2(:)]';
            cell_vec = num2cell(vec,[1,3]);
            if iscell( format )
                h = plot3( cell_vec{:}, format{:} );
            else
                h = plot3( cell_vec{:}, format );
            end
        end
    end
    
end

