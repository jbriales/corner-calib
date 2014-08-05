function imgtrack = initialisation_calib(img_gray)

% Initialisation of the corner tracking algorithm

% Input:
% img_gray  - Gray image in the first iteration

% Outputs:
% x0        - Initial estimation [px py theta_1 theta_2 theta_3]
% Q         - Initial end points (for X, Y and Z axes)

% Grab the points of the Y-junction (first the corner)
h = figure; imshow(img_gray);
title('Grab the center, X, Y and Z points')
[x, y] = getpts;
close(h);

% End points assignment
Q = [ x(2:end), y(2:end) ]';

% Estimate the initial values
px = x(1); py = y(1); 
p_0 = [px; py];
for i = 1:3
%     N_k(i,1) = x(i+1)-xmin +1;
%     N_k(i,2) = y(i+1)-ymin +1;
    v(:,i) = Q(:,i) - p_0;
    v(:,i) = v(:,i)./norm(v(:,i));    
%     n_k(i,1) = -v(i,2);
%     n_k(i,2) = v(i,1);
end
t1 = atan2(v(2,1),v(1,1));
t2 = atan2(v(2,2),v(1,2));
t3 = atan2(v(2,3),v(1,3));

x0 = [px py t1 t2 t3]';     % Initial solution

imgtrack = struct('x', x0,...
                  'mag_th', [0 0 0],...
                  'ang_th', [0 0 0]);
end