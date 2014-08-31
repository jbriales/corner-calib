function scans = loadLidarData( typeOfSource, path, tag )
% scan = readLidarData( typeOfSource, path )
%   typeOfSource - format of stored data (could be Rawlog or Blender)

if ~exist('path','var')
    path = uigetdir('','Choose dataset parent folder');
end

typeOfSource = lower( typeOfSource );
switch typeOfSource
    case 'rawlog'
        scans = loadRawlog( fullfile( path ), tag );
    case 'blender'
        scans = loadBlender( fullfile( path ) );
    otherwise
        error('Type of source specified does not exist')
end

end

function scans = loadRawlog( path, tag )
% scans = loadRawlog( path, tag )
% Read data from rawlog dataset (scan)
% Input:
%   path - path to dataset parent folder 
%   tag  - laser files base tag (after rawlog)
% Output:
%   scans - Array with structs containing:
%       xy   - 2xN array with all detected points
%       mask - validity (logical) mask (0 for wrong points)
%       ts   - vector with timesample values
%   folder - path to selected dataset

% Load range data and times
% TODO: Change rawlog obtention
% error( 'Adapt rawlog value to be more flexible' )
[~,rawlog]  = fileparts( path );
% tag = '_LASER_HOKUYO_UTM';
file_range  = fullfile( path, 'laser', strcat(rawlog, tag, '.txt') );
file_times  = fullfile( path, 'laser', strcat(rawlog, tag, '_times.txt') );
file_config = fullfile( path, 'laser', strcat(rawlog, tag, '_config.txt') );
R  = load( file_range );
ts = load( file_times );

Nobs = size(ts,1);

scan       = readConfigFile( file_config );
scan.FOV   = deg2rad(scan.FOVd);
scan.theta = linspace( -scan.FOV/2, +scan.FOV/2, scan.Npts );
scan.cth   = cos( scan.theta );
scan.sth   = sin( scan.theta );

x = repmat(scan.cth, Nobs,1) .* R;
y = repmat(scan.sth, Nobs,1) .* R;

% Remove zero entries (mask 1 for correct values)
mask_zero = R > 0;
% Crop scan with given angle (mask 1 for values inside)
mask_crop = abs(scan.theta) < deg2rad( scan.crop / 2 );

% Output mask (mask 1 for valid values (columns) in x_scan and y_scan)
mask = mask_zero( :, mask_crop );

x = x( :, mask_crop );
y = y( :, mask_crop );

Nsam = size(mask_crop,2); % TODO: CHECK value
voidStruct = struct('xy',  zeros(2,Nsam),...
                    'mask',zeros(1,Nsam),...
                    'ts',  [] );
scans = repmat( voidStruct, Nobs, 1 );

for i=1:Nobs
    scans(i).xy   = [ x(i,:) ; y(i,:) ];
    scans(i).mask = mask(i,:);
    scans(i).ts   = ts(i);
    scans(i).metafile = fullfile(path,'meta_laser',num2str(ts(i),'%.6f'));
end

end
       

function scans = loadBlender( path )

files = dir( fullfile( path, 'laser' ) );
files(1:2) = []; % Deletes first two elements (current and parent directory)
files = {files(:).name}; % Keeps only the names of files

Nobs = length(files);
voidStruct = struct('xy',  [],...
                    'x', [],...
                    'y', [],...
                    'mask',[],...
                    'ts',  [] );
scans = repmat( voidStruct, Nobs, 1 );
for i=1:Nobs
    name = files{i};
    
    ts = sscanf(name,'%f.pcd');
    pts = double(loadpcd( fullfile(path,'laser',name) ));
    
    % Scan frame needs to be transformed from Blender to typical
    % Scanner-to-World coordinates
    R_s_b = [ 0 0 -1
             -1 0  0
             0 1  0 ];
    pts = R_s_b * pts(1:3,:);
    
%     scans(i).x  = pts(1,:);
%     scans(i).y  = pts(2,:);
    scans(i).xy   = pts(1:2,:);
    scans(i).mask = ones(1, size(pts,2));
    scans(i).ts   = ts;
    scans(i).metafile = num2str(ts,'%.6f');
    scans(i).metafile = fullfile(path,'meta_laser',num2str(ts(i),'%.6f'));
end

end