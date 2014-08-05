%% Simulator for the analysis of the errors in noisy cases
% -------------------------------------------------------

% Previous assertions:
clear; clc;
close all;
%dbstop if error;
kbhit('stop');
kbhit('init');      % For stopping anywhere after keyboard hit

% Control variables
mask = [0 0.001 0.01 0.1 1.0];
N    = 20:10:100;               % number of frames for each simulation
K    = 100;                       % number of simulations

% Simulation parameters
sd_pixel    = 0.5;    % standard deviation of the camera
sd_laser    = 0.01;   % standard deviation of the laser
outlier_k   = 0;      % proportion of outliers for each simulation

%% FOR LOOP BEGIN
% Simulate for different noise in the images and laser, for a variable
% number of frames
sd_pixel_vec = mask * sd_pixel;
sd_laser_vec = mask * sd_laser;

for n = 1:size(N,2)
    m = N(n);
    for i = 1 : size(mask,2)
        sd_pixel = sd_pixel_vec(i);
        for j = 1 : size(mask,2)
            sd_laser = sd_laser_vec(j);   
            for k = 1 : K
                % Simulation saving the data in the structure 'results'
                results(n,i,j,k) = simulator( m, sd_pixel, sd_laser, outlier_k );    
            end           
        end
    end
end

save('simulation.m');

% %% PLOT THE RESULTS
% % TODO: rest of the variables expressed in the title of the graphic
% i = 2;
% j = 2;
% plot_N(results, N, i, j );
% plot_pixel(results, sd_pixel_vec, i, j );
% plot_laser(results, sd_laser_vec, i, j );