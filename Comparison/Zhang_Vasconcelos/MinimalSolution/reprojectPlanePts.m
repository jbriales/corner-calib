% x = [r t]; xdata = [Kvec P]

function y = reprojectPlanePts(x,xdata)

% get homography from plane pose
R = aa2R([x(1);x(2);x(3)]);
H = [R(1:3,1:2) [x(4);x(5);x(6)]];   

% % Option 1: Tune camera intrinsics
% f  = x(7);
% ar = x(8);
% sk = x(9);
% c1 = x(10);
% c2 = x(11);
% K  = [f*ar sk*f c1; 0 f/ar c2; 0 0 1];
% P  = reshape(xdata,2,[]);
% P  = [P; ones(1,size(P,2))];

% Option 2: Do not tune camera intrinsics
f  = xdata(1);
ar = xdata(2);
sk = xdata(3);
c1 = xdata(4);
c2 = xdata(5);
K  = [f*ar sk*f c1; 0 f/ar c2; 0 0 1];

P = reshape(xdata(6:end),2,[]);
P = [P; ones(1,size(P,2))];

y = H*P;

y = y./(ones(3,1)*y(3,:));
y = K * y;% AddDistDivModel(y, eye(3), -0.3994);

y = eye(2,3) * (y./(ones(3,1)*y(3,:)));



