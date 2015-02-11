% Script to test trel

% Add necessary toolboxes
addpath('/home/jesus/Dropbox/Computer_Vision/Matlab_tests/VisionBase');
addpath(genpath('/home/jesus/Dropbox/Computer_Vision/Matlab_tests/PeterKoversi'));
cd /home/jesus/Code/Corner_Calibration
addpath(genpath(pwd))

test = COdometryTest(false); test.WITH_DEBUG = false;
test.Cam.sd = 0.1;

arr_x = logspace(-2,0,15);
Nsim = 1e3;

test.test( 'trel_sd', arr_x, Nsim );
