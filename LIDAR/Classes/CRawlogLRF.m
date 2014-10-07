classdef CRawlogLRF
    %CRawlogLRF Container for Rawlog laser data
    %   The rawlog LRF data should be previously extracted through
    %   RawlogViewer into txt files
    % The files structure is:
    % dataset_date_tag(_file).txt
    %   dataset - dataset name (typically dataset)
    %   date    - date and time of dataset
    %   tag     - sensorLabel given in .ini
    %   file    - could be empty (for ranges), _times
    
    properties
        scans % Array of scans
        
        % Useful data
        Nobs
        
        % Auxiliary rawlog data
        dataset
        date
        tag
    end
    
    methods
        function obj = CRawlogLRF( folder, config )
            % Get directory files
            files = dir(folder);
            files(1:2) = [];
            
            % Extract filename parts of range file (basename) and use to take the rest:
            [~,name,ext] = fileparts( files(1).name );
            file_range = fullfile( folder, strcat( name, ''      , ext ) );
            file_times = fullfile( folder, strcat( name, '_times', ext ) );
            
            rr = load( file_range );
            ts = load( file_times );
            
            Nobs = numel( ts );
            
            % Store object properties
            obj.Nobs = Nobs;
            obj.dataset = []; % Read from auxiliary file?
            obj.date = [];
            obj.tag = [];
            
            % Create array of scans objects and store (in parallel
            % assignment with deal)
            cell_rr = num2cell( rr, 2 );
            cell_ts = num2cell( ts );
            obj.scans = CScan( cell_rr, cell_ts );
        end
        
    end
    
end

