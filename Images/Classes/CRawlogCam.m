classdef CRawlogCam
    %CRawlogCam Class for interaction with stored image datasets
    %   This class aims to provide seamless access to stored image data
    %   isolating from dataset type
    
    properties
        frames % Cell array of images
        
        % Useful data
        Nobs
        
        % Auxiliary rawlog data
        typeOfSource
        dataset
        date
        tag
        
        % Other options
        preload = false
    end
    
    methods
        function this = CRawlogCam( folder, config, typeOfSource, preload )
            % this = CRawlogCam( folder, config, typeOfSource, preload )
            % Constructor for CRawlogCam
            % Input:
            % FOLDER is the base container directory for the dataset files
            % and subfolders
            % CONFIG is the configuration object for the camera
            % TYPEOFSOURCE is the format of dataset directory
            % PRELOAD selects if preloading all images            
            
            switch typeOfSource
                case 'Freiburg'
                    rgb_dir = fullfile( folder, 'rgb' );
                    files = dir(rgb_dir);
                    files(1:2) = [];
                    
                    this.Nobs = numel( files );
                    % Initialize frames
                    
                    this.frames = CFrame.empty( 0, this.Nobs );
                    for i=1:numel( files )
                        this.frames(i) = CFrame( ...
                            fullfile( rgb_dir, files(i).name ) );
                    end
                    
                    % Dataset info
                    this.typeOfSource = typeOfSource;
                    % Get dataset folder name
                    [pathstr,name,~] = fileparts(folder);
                    if isempty(name)
                        % Repeat if name empty (case of ending /)
                        [pathstr,name,~] = fileparts(pathstr);
                    end
                    this.dataset = name;
            end
        end
    end
    
end

