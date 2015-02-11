% Script to test perspective distortion

% Add necessary toolboxes
addpath('/home/jesus/Dropbox/Computer_Vision/Matlab_tests/VisionBase');
addpath(genpath('/home/jesus/Dropbox/Computer_Vision/Matlab_tests/PeterKoversi'));
cd /home/jesus/Code/Corner_Calibration
addpath(genpath(pwd))

test = COdometryTest(false); test.WITH_DEBUG = false;
test.Cam.sd = 0.1;

arr_x = logspace(-3,log10(0.5),20);
Nsim = 1e2;

test.test( 'persp', arr_x, Nsim );
