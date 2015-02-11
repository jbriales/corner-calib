classdef CFrame < handle
    %CFrame Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Image data
        img
        
        % Pose data (if available)
        pose
        
        % Timesample
        ts
        
        % Storage information
        path
        metapath
        file
        metafile
        file_idx
        
        preload = false
    end
    
    properties (Dependent)
        path_file
        path_metafile
    end
    
    methods
        function this = CFrame( S )
            if nargin~=0
                if isstruct( S )
                    % Array constructor from structure
                    this(1,numel(S)) = CFrame;
                    C = fieldnames( S );
                    C = { C{:} };
                    for field = C
                        field = field{:}; %#ok<FXSET>
                        if isprop(this,field)
                            [this.(field)] = deal(S.(field));
                        end
                    end
                elseif ischar( S )
                    % Constructor from file path
                    path_and_file = S;
                    [pathstr,name,ext] = fileparts(path_and_file);
                    
                    this.ts = str2double( name );
                    this.path = pathstr;
                    this.metapath = [pathstr,'-meta'];
                    this.file = [name,ext];
                    this.metafile = [name,ext];
                else
                    error('Not valid input for this class');
                end
                % Assure existence of metadata folder
                if ~exist(this.metapath,'dir')
                    mkdir( this.metapath );
                end
            end
        end
        
        function path_file = get.path_file( this )
            path_file = fullfile( this.path, this.file );
        end
        function path_metafile = get.path_metafile( this )
            path_metafile = fullfile( this.metapath, this.metafile );
        end
        
        function img = loadImg( this )
            % WARNING-TODO: Previous use included img, change path input to
            % rgb or img folder
%             img = imread( fullfile(this.path,'img',this.file) );
            img = imread( fullfile(this.path,this.file) );
        end
        function deleteImg( this )
            file = fullfile(this.path,'img',this.file);
            if exist( file, 'file' )
                answer = input('The image will be deleted: (y/n)','s');
                if strcmpi(answer,'Y')
                    delete( file );
                    warning('Image %s removed',file);
                end
            end
        end
        
    end
    
end

