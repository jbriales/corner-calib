function imgs = loadCamData( path, imgFormat, ptsFormat, blur, decimation )
% imgs = loadCamData( path, imgFormat, ptsFormat )
%   path - path of the parent dataset directory
%   imgFormat - string with format of input images
%   ptsFormat - (TODO) string with point storage files

if ~exist('path','var')
    path = uigetdir('','Choose dataset parent folder');
end

if ~exist('imgFormat','var')
    imgFormat = '_%f.jpg';
end
if ~exist('ptsFormat','var')
    ptsFormat = '_%f.pts';
end
if ~exist('decimation','var')
    decimation = 1;
end

% Choose and read images
K          = load( fullfile( path, 'intrinsic_matrix.txt' ) );
distortion = load( fullfile( path, 'distortion.txt' ) );
[~,~,img_ext] = fileparts(imgFormat);
img_files  = uigetfile( fullfile(path,'img',strcat('*',img_ext)), 'Choose the images to calibrate', 'MultiSelect', 'on' );
if isnumeric( img_files ) && img_files == 0
    imgs = [];
    return
end
img_idxs = 1:decimation:length(img_files);
imgs = repmat( struct('I',[], 'file',[], 'ts',[], 'K',[], 'P',[] ), length(img_idxs), 1 );
for k=1:length(imgs)
    i = img_idxs(k);
    
    Im_ = imread( fullfile(path,'img',img_files{i}) );
    if any( distortion )
        Im_ = cv.undistort(Im_, K, distortion);
    end
    if blur
        G = fspecial('gaussian',[5 5],2);
        %# Filter it
        Im_ = imfilter(Im_,G,'same');
    end
    
    imgs(k).I = Im_;
    imgs(k).file = img_files{i};
    imgs(k).ts = sscanf(img_files{i}, imgFormat);
    imgs(k).K = K;
%     if exist(fullfile(path,'img',strcat(img_files{i}(1:end-3),'pts')),'file')
%         imgs(i).P = load(fullfile(path,'img',strcat(img_files{i}(1:end-3),'pts')));
%     end
end