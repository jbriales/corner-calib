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

LRF.tracking = true;
% LRF.tracking = true;
Cam.tracking = true;
% Cam.tracking = true;
if ~exist('ptime','var')
    ptime = 0.01;
    warning('ptime set to default value: %f',ptime);
end
for nobs=1:length(scans)
    fprintf('Frame %d\n',nobs);
    
    Cam.setFrame( frames(nobs) );
    obj_xi = Cam.computeXi;
    Cam.visualize;
    
    LRF.setFrame( scans(nobs) );
    [v,A_v,q,A_q] = LRF.compute;
    LRF.visualize;
    pause(ptime)
    if upper(kbhit) == 'L' % To stop after pressing 'L' key
        % To fix a wrong scanner segmentation:
        LRF.fixFrame( scans(nobs) );
    end
end