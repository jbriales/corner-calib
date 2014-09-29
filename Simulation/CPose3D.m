classdef CPose3D < handle
    %CPose Pose of a 3D object wrt reference (typically World)
    %   Fully defined through its orientation and 
    %   translation wrt reference frame
    %   Constructor:
    %   pose = CPose3D( R, t )

    properties
        R   % 3x3 rotation matrix of Polygon seen from World
        t   % 3x1 translation vector of Polygon seen from World
    end
    
    properties (SetAccess = protected, Dependent) % (Read-only)
        T       % 4x4 pose of Polygon seen from World
    end
    
    methods
        % Constructor
        function obj = CPose3D( R, t )
            obj.R  = R;
            obj.t  = t;
        end
        
        % Get methods
        function T = get.T( obj )
            T = [ obj.R obj.t ; zeros(1,3) 1 ];
        end
        
        % Plotting methods
        % Plot object frame in 3D space
        function h = plotFrame( obj, label, color ) 
            h = plotframe( obj.T, 1, label, color );
        end
        % Plot reference frame in 3D space (World)
        % Not working?
        function h = plotReference( obj, label, color )
            if ~exist('label','var')
                label = 'W';
            end
            if ~exist('color','var')
                color = 'k';
            end
            h = plotframe( eye(4), 1, label, color );
        end
    end
end