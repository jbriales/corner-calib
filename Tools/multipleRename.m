function multipleRename(STARTPATH)
% Function that rename multiple files in a folder to other folder according
% to a pattern (usually number-container filenames)

% First version: To be improved

if ~exist('STARTPATH','var')
    STARTPATH = '';
end

% Get root directory where two folders (old and new) are contained
input_fold = uigetdir(STARTPATH,'Choose source directory'); % CHANGE THIS VALUE
output_fold = uigetdir(STARTPATH,'Choose destination directory'); % CHANGE THIS VALUE

if isnumeric( input_fold ) || isnumeric( output_fold )
    if input_fold == 0 || output_fold == 0
        warning('Select two (different) folders')
        return
    end
end

im_name = dir( input_fold );
im_name(1:2) = []; % Deletes first two elements
im_name = {im_name(:).name}; % Keeps only the names of im_name

for i=1:length(im_name)
    name1 = im_name{i};
    
    number = sscanf(name1,'%d.txt');
    
    name2 = sprintf('%.3d.txt',number+1);
    
    copyfile( fullfile(input_fold,name1), fullfile(output_fold,name2) )
end