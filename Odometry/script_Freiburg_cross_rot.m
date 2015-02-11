% Script to do test on Freiburg data, taking all possible combinations

% Add necessary toolboxes
addpath('/home/jesus/Dropbox/Computer_Vision/Matlab_tests/VisionBase');
addpath(genpath('/home/jesus/Dropbox/Computer_Vision/Matlab_tests/PeterKoversi'));
cd /home/jesus/Code/Corner_Calibration
addpath(genpath(pwd))

test = COdometryTest(false); test.WITH_DEBUG = false;

test.testOrientation( '/media/storage/Research/Dataset/rgbd_dataset_freiburg3_cabinet_final',...
        'Freiburg', [1 60], 30, false );