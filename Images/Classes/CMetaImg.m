classdef CMetaImg < handle
%CMetaImg Summary of this class goes here
% Detailed explanation goes here
properties
xi % Xi object representing corner structure
c % Vertex point
q % Extreme points
r % Length of segments

v_c
v_ang
% Image analysis data
mag_th
ang_th
% Control parameters
optimized % Boolean 1 if optimization has been already applied to xi
end
methods
function obj = CMetaImg( c, q, mag_th, ang_th )
if nargin ~= 0
obj.c = c;
obj.q = q;
[v,r] = deal(cell(1,3));
for k=1:3
v{k} = snormalize( q{k} - c );
r{k} = norm( q{k} - c );
end
obj.r = r;
obj.xi = Cxi( c, v{1}, v{2}, v{3} );
obj.mag_th = mag_th;
obj.ang_th = ang_th;
end
end
function manualHint( obj )
% manualHint( obj )
% Get image data for tracking by manually clicking points
% It should be applied after showing image in some figure
% Outputs:
% x0 - Initial estimation [px py theta_1 theta_2 theta_3]
% Q - Initial end points (for X, Y and Z axes)
title(sprintf('Grab the center, X, Y and Z points\tIMPORTANT: Make coordinates dextrorotatory!'))
[x, y] = getpts;
while isempty(x)
warning('Input points');
[x,y] = getpts;
end
% Vertex point
c = [ x(1), y(1) ]';
% End points assignment
Q = [ x(2:end), y(2:end) ]';
Q = num2cell( Q, 1 );
for k = 1:3
v{k} = snormalize( Q{k} - c );
r{k} = norm( Q{k} - c );
end
% Create xi object from above data
xi = Cxi( c, v{1}, v{2}, v{3} );
% Set object properties values
obj.c = c;
obj.q = Q;
obj.r = r;
obj.xi = xi;
% Set optimized property as false (it is necessary to do refinement)
obj.optimized = false;
end
end
end
