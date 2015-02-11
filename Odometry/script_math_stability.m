% Script to test numerical stability of P3oA solution wrt Ali

% Add necessary toolboxes
addpath('/home/jesus/Dropbox/Computer_Vision/Matlab_tests/VisionBase');
addpath(genpath('/home/jesus/Dropbox/Computer_Vision/Matlab_tests/PeterKoversi'));
cd /home/jesus/Code/Corner_Calibration
addpath(genpath(pwd))

test = COdometryTest(false); test.WITH_DEBUG = false;
test.Cam.sd = 0;

Nsim = 1e5;

test.test_Math( Nsim )
