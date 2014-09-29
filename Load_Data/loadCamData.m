function imgs = loadCamData( path, imgFormat, blur, decimation, idx_span )
% imgs = loadCamData( path, imgFormat, ptsFormat )
%   path - path of the parent dataset directory
%   imgFormat - string with format of input images

if ~exist('path','var')
    path = uigetdir('','Choose dataset parent folder');
end

if ~exist('imgFormat','var')
    imgFormat = '_%f.jpg';
end
if ~exist('decimation','var')
    decimation = 1;
end

% Choose and read images
K          = load( fullfile( path, 'intrinsic_matrix.txt' ) );
distortion = load( fullfile( path, 'distortion.txt' ) );
[~,~,img_ext] = fileparts(imgFormat);

if ~exist('idx_span','var') || isempty(idx_span)
    img_files = uigetfile( fullfile(path,'img',strcat('*',img_ext)), 'Choose the images to calibrate', 'MultiSelect', 'on' );
    if isnumeric( img_files ) && img_files == 0
        imgs = [];
        return
    end
    img_idxs = 1:decimation:length(img_files);
else
    img_files = dir(fullfile(path,'img'));
    img_files(1:2) = []; % . and .. items
    img_files = {img_files.name};
    % Filter with label
    prefix = imgFormat( 1:strfind(imgFormat,'%')-1 );
    for k=length(img_files):-1:1
        if isempty( strfind( img_files{k}, prefix ) )
            img_files(k) = [];
        end
    end
    img_idxs = idx_span(1):decimation:idx_span(end);
end
imgs = repmat( struct('I',[], 'file',[], 'path',[], 'ts',[], 'K',[], 'distortion', [] ), length(img_idxs), 1 );
for k=1:length(imgs)
    i = img_idxs(k);
    
    [~,img_name,~] = fileparts( img_files{i} );
    metafile = strcat( img_name, '.mat' );
    
    % Deprecated: Images are loaded only at time of use
%     Im_ = imread( fullfile(path,'img',img_files{i}) );
%     if any( distortion )
%         Im_ = cv.undistort(Im_, K, distortion);
%     end
%     if blur
%         G = fspecial('gaussian',[5 5],2);
%         %# Filter it
%         Im_ = imfilter(Im_,G,'same');
%     end
    
%     imgs(k).I = Im_;
    imgs(k).file = img_files{i};
    imgs(k).metafile = metafile;
    imgs(k).file_idx = int8(img_idxs(k));
%     imgs(k).path = fullfile(path,'img',img_files{i});
    imgs(k).path = fullfile(path);
%     imgs(k).metafile = fullfile(path,'meta_img',img_files{i});
    imgs(k).ts = sscanf(img_files{i}, imgFormat);
    imgs(k).K = K;
    imgs(k).distortion = distortion;
end