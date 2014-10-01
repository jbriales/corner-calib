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

Cam.tracking = true;
% Cam.tracking = true;
Cam.linesToEnd = linesToEnd;
% Cam.linesToEnd = [false, true, true];
% Cam.linesToEnd = [true, true, true];
% Cam.linesToEnd(1) = ~Cam.linesToEnd(1);
% Cam.linesToEnd(2) = ~Cam.linesToEnd(2);
% Cam.linesToEnd(3) = ~Cam.linesToEnd(3);
if ~exist('ptime','var')
    ptime = 0.1;
    warning('ptime set to default value: %f',ptime);
end
for nobs=1:length(scans)
    try
        fprintf('Frame %d\n',img_idxs(1)+nobs);
        Cam.setFrame( frames(nobs) );
        %     Cam.visualize; % Current image with previous lines
        obj_xi = Cam.computeXi;
        Cam.visualize;
        % Optional for debug:
        % hInlPts = Cam.GOptimizer.plotInlierPts
        % hInlW   = Cam.GOptimizer.plotWeightPts
        pause(ptime)
        if upper(kbhit) == 'L' % To stop after pressing 'L' key
            Cam.fixFrame( frames(nobs) );
        end
%         kh = upper(kbhit);
%         if ~isempty(kh)
%             switch kh
%                 case 'D'
%                     ptime = ptime * 2;
%                 case 'A'
%                     ptime = ptime / 2;
%                 case 'S'
%                     Cam.fixFrame( frames(nobs) );
%             end
%             clc
%             figure(gcf)
%         end
    % frames(nobs).deleteImg;
    catch exception
        warning(exception.message)
        keyboard
        Cam.fixFrame( frames(nobs) );
    end
end