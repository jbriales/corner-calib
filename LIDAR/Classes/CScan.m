classdef CScan < handle
    %CScan Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Minimal data (as stored in Rawlog)
        r   % Range data (in m)
        ts  % timesample
        
        % Scan data
        xy
        mask
        
        % Some extra options
        crop
        
        % Timesample
        delta_ts % Difference with associated data (camera)
        
        % Storage information
        path
        metafile
    end
    
    methods
        function obj = CScan( in1, ts )
            % obj = CScan( S )
            % obj = CScan( r, ts )
            if nargin~=0
                switch nargin
                    case 1
                        S = in1;
                    case 2
                        r = in1; %#ok<*PROP>
                end
                if exist('r','var')
                    % Array constructor
                    if iscell( r )
                        obj(1,numel(ts)) = CScan;
                        [obj.r] = deal( r{:} );
                        [obj.ts] = deal( ts{:} );
                    else
                        obj.r = r;
                        obj.ts = ts;
                    end
                    
                elseif exist( 'S', 'var' )
                    % Array constructor
                    obj(1,numel(S)) = CScan;
                    C = fieldnames( S );
                    C = { C{:} };
                    for field = C
                        field = field{:}; %#ok<FXSET>
                        if isprop(obj,field)
                            [obj.(field)] = deal(S.(field));
                        end
                    end
                else
                    error('Not valid input for this class');
                end
            end
        end
        
        function loadScan( obj )
            % Blender case
            obj.xy = loadBlenderPCD( obj.path );
        end 
    end
    
end

