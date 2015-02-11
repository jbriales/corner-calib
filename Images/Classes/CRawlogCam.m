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
                    %% FREIBURG
                    % First check if the dataset was already processed
                    if exist(fullfile(folder,'dataset.mat'),'file')
%                         yn = input('Use previous loaded data? ');
                        fprintf('Using previous loaded data in file:\n%s\n',...
                                fullfile(folder,'dataset.mat'));
                        load( fullfile(folder,'dataset.mat') );
                        return
                    end
                    
                    rgb_dir = fullfile( folder, 'rgb' );
                    files = dir(rgb_dir);
                    files(1:2) = [];
                    
                    this.Nobs = numel( files );
                    % Read image file data (without loading image)
                    this.frames = CFrame.empty( 0, this.Nobs );
                    for i=1:numel( files )
                        this.frames(i) = CFrame( ...
                            fullfile( rgb_dir, files(i).name ) );
                    end
                    
                    % Get GT pose from Freiburg dataset
                    FID = fopen(fullfile(folder,'groundtruth.txt'));
                    fgetl( FID ); fgetl( FID ); fgetl( FID ); % Read the first 3 rows
                    v = fscanf(FID, '%f %f %f %f %f %f %f %f\n', [8 Inf]);
                    fclose(FID);
                    ts = v(1,:); % Timestamps vector
                    % Assign each frame the closest GT value
                    for i=1:numel( files )
                        [~, gt_i] = min( abs( this.frames(i).ts - ts ) );
                        T_ = quaternion2matrix( [v(8,gt_i), v(5:7,gt_i)'] );
                        t_ = v(2:4,gt_i); % Translation [tx,ty,tz]'
                        this.frames(i).pose = CPose3D( T_(1:3,1:3), t_ );
                    end
                    clear gt_i T_ t_ % Dispose of temporal variables
                    
                    % Dataset info
                    this.typeOfSource = typeOfSource;
                    % Get dataset folder name
                    [pathstr,name,~] = fileparts(folder);
                    if isempty(name)
                        % Repeat if name empty (case of ending /)
                        [pathstr,name,~] = fileparts(pathstr);
                    end
                    this.dataset = name;
                    
                    % Save dataset for further uses
                    save( fullfile(folder,'dataset.mat'), 'this' );
                    
                case 'Blender'
                    %% BLENDER
                    % First check if the dataset was already processed
                    if exist(fullfile(folder,'dataset.mat'),'file')
%                         yn = input('Use previous loaded data? ');
                        warning('Using previous loaded data in file:\n%s\n',...
                                fullfile(folder,'dataset.mat'));
                        load( fullfile(folder,'dataset.mat') );
                        return
                    end
                    
                    rgb_dir = fullfile( folder, 'rgb' );
                    files = dir(rgb_dir);
                    files(1:2) = [];
                    
                    this.Nobs = numel( files );
                    % Read image file data (without loading image)
                    this.frames = CFrame.empty( 0, this.Nobs );
                    for i=1:numel( files )
                        this.frames(i) = CFrame( ...
                            fullfile( rgb_dir, files(i).name ) );
                    end
                    
                    % Get GT pose from Freiburg dataset
                    FID = fopen(fullfile(folder,'groundtruth.txt'));
                    fgetl( FID ); fgetl( FID ); fgetl( FID ); % Read the first 3 rows
                    v = fscanf(FID, '%f %f %f %f %f %f %f %f\n', [7 Inf]);
                    fclose(FID);
                    ts = v(1,:); % Timestamps vector
                    % Assign each frame the closest GT value
                    for i=1:numel( files )
                        [ts_err, gt_i] = min( abs( this.frames(i).ts - ts ) );
                        
                        r_xyz = v(5:7,gt_i);
                        % IMPORTANT:
                        % Change sign of Y and Z rotations when correction is needed from Blender convention
                        % In Blender the camera is pointing in the opposite direction to Z axis
                        R_ = RotationZ( deg2rad(r_xyz(3)) ) *...
                              RotationY( deg2rad(r_xyz(2)) ) *...
                              RotationX( deg2rad(r_xyz(1)) ) *...
                              RotationX( pi ); % 180[deg] wrt X axis
%                         r_xyz(2:3) = -r_xyz(2:3);
%                         R_ = RotationYPR( deg2rad( r_xyz(end:-1:1) ) ) ;
                        t_ = v(2:4,gt_i); % Translation [tx,ty,tz]'
                        this.frames(i).pose = CPose3D( R_, t_ );
                    end
                    clear gt_i R_ t_ % Dispose of temporal variables
                    
                    % Dataset info
                    this.typeOfSource = typeOfSource;
                    % Get dataset folder name
                    [pathstr,name,~] = fileparts(folder);
                    if isempty(name)
                        % Repeat if name empty (case of ending /)
                        [pathstr,name,~] = fileparts(pathstr);
                    end
                    this.dataset = name;
                    
                    % Save dataset for further uses
                    save( fullfile(folder,'dataset.mat'), 'this' );
            end
        end
    end
    
end

