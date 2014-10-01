classdef CScan < handle
    %CScan Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Scan data
        xy
        mask
        
        % Timesample
        ts
        delta_ts
        
        % Storage information
        path
        metafile
    end
    
    methods
        function obj = CScan( S )
            if nargin~=0
                if isstruct( S )
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

