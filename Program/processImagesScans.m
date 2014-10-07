% Previous assertions:
clear
close all
dbstop if error % Set debug mode on if error occurs
kbhit('stop')
kbhit('init') % Track keyboard hits

% To position figures centered and visible by default
set(0,'defaultfigureposition', [642   503   560   420])

main_type = readConfigFile( fullfile(pwd,'main.ini'),'','main_type' );
userOpts = readConfigFile( fullfile(pwd,'main.ini'), main_type );
extractStructFields( userOpts );

all_dataset_data = loadDataset( userOpts );
extractStructFields( all_dataset_data );

set(gcf,'Position', videoFigurePos);
% tightfig
videofile = fullfile(path,'Video.avi');
% videofile = fullfile(path,'Video');
if exist(videofile,'file')
    video = VideoReader(videofile);
end
nvideoframe = 330;
% im = read(video,k);
t0 = frames(1).ts;

video_writer = VideoWriter('video_out.avi');
open(video_writer);

LRF.tracking = true;
% LRF.tracking = true;
Cam.tracking = true;
% Cam.tracking = true;
if ~exist('ptime','var')
    ptime = 0.01;
    warning('ptime set to default value: %f',ptime);
end
for nobs=1:length(scans)
    tic;
    fprintf('Frame %d\n',nobs);
    
    Cam.setFrame( frames(nobs) );
    obj_xi = Cam.computeXi;
    Cam.visualize;
    title('Camera capture','FontSize',20);
    
    LRF.setFrame( scans(nobs) );
    [v,A_v,q,A_q] = LRF.compute;
    LRF.visualize;
%     LRF.visualize( LRF_video_borders );
    title('LRF 2D scan','FontSize',20);
    % Wait to frequency time
    if 0
        t1 = toc;
        while( t1 < 1/30 )
            t1 = toc;
            disp(toc);
        end
        disp(1/toc);
    else
        pause(0.01)
    end
 
    % Show video
    inc_frame = (frames(nobs).ts - t0)*30;
    videoimg = read(video,nvideoframe+inc_frame);
    subplot(hVid);
    imshow(videoimg);
%     nvideoframe = nvideoframe + 1;
    
if 1
    file_out = fullfile(path,'VideoOut',Cam.frame.file);
    hgexport(gcf, file_out, hgexport('factorystyle'), 'Format', 'jpeg');
end
    
    if upper(kbhit) == 'L' % To stop after pressing 'L' key
        % To fix a wrong scanner segmentation:
        LRF.fixFrame( scans(nobs) );
    end
end