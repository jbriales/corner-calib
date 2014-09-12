% x = [rlrs2cam tlrs2cam rp2cam1 tp2cam1 ... rp2camN tp2camN]
% xdata = [Kvec nPlanePts nLrsPts planePts lrsPts]
% y = [camPts1 ... camPtsN 0 0 ... 0]

function y = reprojectLaserCamCalib(x,xdata)

rlrs2cam = [x(1); x(2); x(3)];
tlrs2cam = [x(4); x(5); x(6)];
Rlrs2cam = aa2R(rlrs2cam);
Tlrs2cam = [Rlrs2cam tlrs2cam; 0 0 0 1];
p2cam   = x(7:end);
nPlanes = length(p2cam)/6;
if nPlanes ~= floor(nPlanes)
    error('Wrong number of elements in x (should be 6+6*n, n=number of planes)');
end

f  = xdata(1);
ar = xdata(2);
sk = xdata(3);
c1 = xdata(4);
c2 = xdata(5);
Kvec = [f ar sk c1 c2];

nPlanePts = xdata(6:6+nPlanes-1);
nLrsPts   = xdata(6+nPlanes:6+2*nPlanes-1);

planePtsIdx = 6+2*nPlanes;
lrsPtsIdx   = planePtsIdx + 2*sum(nPlanePts);

ycam = [];
ylrs = [];
for i=1:nPlanes  
    clear L
    
    rp2cam = [x(6*i+1); x(6*i+2); x(6*i+3)];
    tp2cam = [x(6*i+4); x(6*i+5); x(6*i+6)];    
    Rp2cam = aa2R(rp2cam);
    Tp2cam = [Rp2cam tp2cam; 0 0 0 1];
    
    % Reproject plane points in camera view
    PlanePts    = xdata(planePtsIdx:planePtsIdx+2*nPlanePts(i)-1);
    planePtsIdx = planePtsIdx + 2*nPlanePts(i);
    data        = [Kvec PlanePts];
    rt          = [rp2cam tp2cam];
    ycam        = [ycam reshape(reprojectPlanePts(rt,data),1,[])];
    
    % Get laser points
    LrsPts    = reshape(xdata(lrsPtsIdx:lrsPtsIdx+3*nLrsPts(i)-1),3,[]);
    lrsPtsIdx = lrsPtsIdx + 3*nLrsPts(i);      
    
    % Get Plane in Laser Frame
    PIlrs = [Rlrs2cam [0;0;0]; -tlrs2cam.'*Rlrs2cam 1]*[Rp2cam [0;0;0]; -tp2cam.'*Rp2cam 1]*[0; 0; 1; 0];
    % Get Laser Lines
    L(1:3,:) = [LrsPts(1:2,:); zeros(1,size(LrsPts,2))];
    for j=1:size(LrsPts,2)
        L(4:6,j) = cross(LrsPts(:,j),[0;0;10]);
    end
    % Get intersection between lines and plane
    P=IntersectionLinePlane(PIlrs, L);
    % Get error between intersections and laser points
    err = sqrt(sum((LrsPts(1:2,:,1)-P(1:2,:)).^2));
    
    ylrs   = [ylrs err];
end

alpha = 1;

y = [alpha*ycam ylrs];