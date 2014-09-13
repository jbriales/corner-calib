%
% [R, m] = rotation3Pt(x1, x2, t, disp)
%
% determines the rotation between 2 camera frames, given the camera 
% projection of 3 points and the translation between the cameras.
%
% INPUT:
%   x1 - 3x3 matrix with columns as point projections in camera 1
%   x2 - 3x3 matrix with columns as point projections in camera 2
%   t  - 3x1 translation vector of camera 1 in camera 2 frame
%
% OUTPUT:
%   R - 3x3xN matrix with N solutions to rotation mapping points from 
%       camera 1 to camera 2
%   m - 4x1xN vector with N solutions to plane containing the 3 projected 
%       points
%   

function [R, m] = rotation3Pt(x1, x2, t, display)

x1(:,1) = x1(:,1)./norm(x1(:,1));
x1(:,2) = x1(:,2)./norm(x1(:,2));
x1(:,3) = x1(:,3)./norm(x1(:,3));

l12 = cross(x1(:,1), x1(:,2));% x1(1:3,1) - x1(1:3,2));
l13 = cross(x1(:,1), x1(:,3));% x1(1:3,1) - x1(1:3,3));
l23 = cross(x1(:,2), x1(:,3));% x1(1:3,2) - x1(1:3,3));

k12 = cross(x2(:,1),x2(:,2));
k13 = cross(x2(:,1),x2(:,3));
k23 = cross(x2(:,2),x2(:,3));

x = [k12 k13 k23];
P = [l12/(t.'*l12) l13/(t.'*l13) l23/(t.'*l23)];

if display
    fprintf('\n\nP3P POINT POSITIONS ERROR\n');
end;
[T, s1, s2, s3] = p3p(x,P,1,display);

R = T(1:3,1:3,:);
m = T(1:4,4,:);

h1 = x1 ./ (ones(3,1)*x1(3,:));
h2 = x2 ./ (ones(3,1)*x2(3,:));

if display
    fprintf('\n\nROTATION ERROR FROM ESSENTIAL MATRIX CONDITION\n');  
    for i=1:size(R,3)
        e = (h2(:,1).'*(R(:,:,i)*skew_symetric_v(t))*h1(:,1))^2 + ...
            (h2(:,2).'*(R(:,:,i)*skew_symetric_v(t))*h1(:,2))^2 + ...
            (h2(:,3).'*(R(:,:,i)*skew_symetric_v(t))*h1(:,3))^2;
        disp(['solution ' int2str(i) ': ' num2str(e)])   
    end
end;


