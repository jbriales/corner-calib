%
% [T, s1, s2, s3] = p3p(x, P, opt, display)
%
% Solves the perspective three point problem using Grunert's algorithm.
% If only positive depths are considered, up to 4 solutions are found,
% otherwise, there are up to 8 solutions.
%
% Implementation is based on Haralick et al, "Review and Analysis of 
% Solutions of the Three Point Perspective Pose Estimation Problem", 1994
%
% INPUT:
%   x   - 3x3 matrix with columns as camera point projections,
%   P   - 3x3 matrix with columns as point coordinates in world frame.
%   opt - 0 to consider positive depths only, and 1 to consider negative 
%         depths as well.
%
%   note: for positive depths only, the third row of x should always be 
%   positive.
%
% OUTPUT:
%   T  - 4x4xN matrix with N solutions to camera pose in world frame
%   s1 - 1xN vector with N solutions to depth of P(:,1) in camera frame.
%   s2 - 1xN vector with N solutions to depth of P(:,2) in camera frame.
%   s3 - 1xN vector with N solutions to depth of P(:,3) in camera frame.
%

function [T, s1, s2, s3] = p3p(x, P, opt, display)

if ~exist('opt','var')
    opt = 0;
end

if ~exist('display','var')
    display = 0;
end

% normalize x for norm 1
x(:,1)=x(:,1)./norm(x(:,1));
x(:,2)=x(:,2)./norm(x(:,2));
x(:,3)=x(:,3)./norm(x(:,3));

% reparameterization
cos_alpha = x(:,2).'*x(:,3);
cos_beta  = x(:,1).'*x(:,3);
cos_gamma = x(:,1).'*x(:,2);

a = norm(P(:,2) - P(:,3));
b = norm(P(:,1) - P(:,3));
c = norm(P(:,1) - P(:,2));

% normalize a, b, c 
k = mean([a b c]);
a = a/k;
b = b/k;
c = c/k;

apc = (a^2+c^2)/b^2;
amc = (a^2-c^2)/b^2;
bmc = (b^2-c^2)/b^2;
bma = (b^2-a^2)/b^2;

A0 = (1 + amc)^2 - (4*a^2/b^2)*cos_gamma^2;
A1 = 4*(-amc*(1+amc)*cos_beta + (2*a^2/b^2)*cos_gamma^2*cos_beta - (1-apc)*cos_alpha*cos_gamma);
A2 = 2*(amc^2 - 1 + 2*amc^2*cos_beta^2 + 2*bmc*cos_alpha^2 - 4*apc*cos_alpha*cos_beta*cos_gamma + 2*bma*cos_gamma^2);
A3 = 4*(amc*(1-amc)*cos_beta - (1-apc)*cos_alpha*cos_gamma + (2*c^2/b^2)*cos_alpha^2*cos_beta);
A4 = (amc-1)^2 - (4*c^2/b^2)*cos_alpha^2;

% quartic equation roots
v = roots([A4 A3 A2 A1 A0]);
% v = unique(real(v(abs((imag(v))<10e-4) | (atan2(real(v),imag(v))<10e-4))));
v = unique(real(v));

for i=1:length(v);
    u(i) = ((1 - amc)*v(i)^2 + 2*amc*cos_beta*v(i) - 1 - amc) / (2*v(i)*cos_alpha - 2*cos_gamma);

    % depths
    s1(i) = sqrt(b^2/(1+v(i)^2-2*v(i)*cos_beta));
    s2(i) = u(i)*s1(i);
    s3(i) = v(i)*s1(i);

    % denormalize
    s1(i) = s1(i)*k;
    s2(i) = s2(i)*k;
    s3(i) = s3(i)*k; 
    
    if opt
        s1(length(v) + i) = -s1(i);
        s2(length(v) + i) = -s2(i);
        s3(length(v) + i) = -s3(i);    
    end
    
    % camera pose in world frame
    T(:,:,i) = absoluteOrientation(P,...
        (ones(3,1)*[s1(i) s2(i) s3(i)]).*x);
    if opt
        T(:,:,length(v) + i) = absoluteOrientation(P,...
            (ones(3,1)*[s1(length(v) + i) s2(length(v) + i) s3(length(v) + i)]).*x);
    end
        
end

if display
    for i=1:size(T,3)
        error(i) = mean([norm(s1(i)*x(:,1) - T(1:3,1:4,i)*[P(:,1);1]) ...
            norm(s2(i)*x(:,2) - T(1:3,1:4,i)*[P(:,2);1]) ...
            norm(s3(i)*x(:,3) - T(1:3,1:4,i)*[P(:,3);1])]);
        disp(['Error of solution ' int2str(i) ': ' num2str(error(i))]);    
    end
end