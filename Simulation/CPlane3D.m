classdef CPlane3D < CPose3D
    %CPlane3D Plane placed in 3D space
    %   Plane object is fully defined through its orientation and 
    %   translation wrt reference frame
    %   Constructor:
    %   plane = CPlane( R, t )
   
    properties (SetAccess = protected, Dependent) % (Read-only)
        n       % 3x1 plane normal vector
        plane   % 4x1 plane vector (normal and distance)
        M       % 4x3 matrix for conversion from 2D to 3D frame
    end
    
    methods
        % Constructor
        function obj = CPlane3D( R, t )
            obj = obj@CPose3D( R, t );
        end
        
        % Convert points in 3D to plane 2D frame (checking if in)
        function pts2D = transform3Dto2D( obj, pts3D )
            N = size( pts3D, 2 );
            rel = pts3D - repmat( obj.t, 1, N );
            % Check that all points belong to plane
            prod = obj.n' * rel;
            max_prod = max( abs( prod ) );
            if max_prod > 1e-6 % micron distance
                error('[CPolygon::transform2D] Points outside plane. max(n''·(p3D-t)=%e',max_prod);
            end
            pts2D = obj.R(:,1:2)' * rel;
        end
        
        % Convert points in plane frame to 3D
        function pts3D = transform2Dto3D( obj, pts2D )
            pts3D = makeinhomogeneous( obj.M * makehomogeneous( pts2D ) );
        end
        
        % Get methods
        function n = get.n( obj )
            n = obj.R(:,3);
        end
        function plane = get.plane( obj )
            plane = [ obj.n ; -obj.n' * obj.t ];
        end
        function M = get.M( obj )
            M = obj.T(:,[1 2 4]); % Remove Z column
        end
        
    end
end