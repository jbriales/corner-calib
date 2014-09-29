classdef CFrame < handle
    %CFrame Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Timesample
        ts
        
        % Storage information
        path
        file
        metafile
        file_idx
    end
    
    methods
        function obj = CFrame( S )
            if nargin~=0
                if isstruct( S )
                    % Array constructor
                    obj(1,numel(S)) = CFrame;
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
        
        function img = loadImg( obj )
            img = imread( fullfile(obj.path,'img',obj.file) );
        end
        
        
    end
    
end

