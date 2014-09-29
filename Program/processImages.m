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

Cam.tracking = false;
% Cam.tracking = true;
for nobs=1:length(scans)
    Cam.setFrame( frames(nobs) );
    obj_xi = Cam.computeXi;
    Cam.visualize;
    keyboard
end