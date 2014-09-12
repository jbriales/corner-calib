close all; clear; clc;
load('rot_input_ex.mat');
n_cam       = [rot_input.N];
line_LRF    = [rot_input.l];

dist  = 0.1;    % Simulation distance (rad)
inc   = 0.01;    % Increments

w_gt    = rodrigues(R_c_s);    % GT rotation vector

w_sim   = [w_gt(2)-dist:inc:w_gt(2)+dist; 
           w_gt(3)-dist:inc:w_gt(3)+dist;];           % Simulated rotation vectors
     
[w_y,w_z] = meshgrid(w_sim(1,:), w_sim(2,:));
w_x       = ones(size(w_y,1),size(w_y,1)) * w_gt(3);       
       
%E = zeros( size(w_sim(1,:),2),size(w_sim(1,:),2) );

for i = 1:size( w_sim(1,:) , 2 )
    for j = 1:size( w_sim(2,:) , 2 )
        
        w = [w_x(i,j) w_y(i,j) w_z(i,j)]';
        R = rodrigues(w);
%         residual = dot( n_cam, R(1:3,1:2) * line_LRF, 1 )';
%         E_trihedron(i,j) = residual' * residual;

        x = [R [0.15 0 0]'];
        theta = Rig.Lidar.FOVd / Rig.Lidar.N ;
        %[residual_,~] = Fun_sum( corner_corresp, x, K );
        [residual_,~] = FunW_sum( corner_corresp, x, K, theta );
        E_kwak(i,j) = residual_;
        
    end
end
% 
% surf(w_y, w_z, E_trihedron);
% view([90 90]);
% 
% figure;

surf(w_y, w_z, E_kwak);
view([90 90]);

