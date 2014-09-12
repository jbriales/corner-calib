function [x0,xdata,ydata] = laserCamOptData(Points,CamCalib,Tp2cam,Tlrs2cam)

rlrs2cam = vrrotmat2vec(Tlrs2cam(1:3,1:3));
rlrs2cam = rlrs2cam(1:3)*rlrs2cam(4);
tlrs2cam = Tlrs2cam(1:3,4);
lrs2cam  = [rlrs2cam tlrs2cam.'];

p2cam    = [];
planePts = [];
camPts   = [];
lrsPts   = [];
for i=1:size(Tp2cam,3)
    rp2cam = vrrotmat2vec(Tp2cam(1:3,1:3,i));
    rp2cam = rp2cam(1:3)*rp2cam(4);
    tp2cam = Tp2cam(1:3,4,i);
    p2cam  = [p2cam rp2cam tp2cam.'];
    
    planePts   = [planePts reshape(Points(i).plane(1:2,:),1,[])];
    camPts     = [camPts reshape(Points(i).camNoisy(1:2,:),1,[])];
    nCamPts(i) = size(Points(i).camNoisy,2);
    lrsPts     = [lrsPts reshape([Points(i).lrsNoisy(1:2,:); zeros(1,size(Points(i).lrsNoisy,2))],1,[])];
    nLrsPts(i) = size(Points(i).lrsNoisy,2);
end

Kvec = [CamCalib.focal CamCalib.aratio CamCalib.skew CamCalib.center(1) CamCalib.center(2)];

x0    = [lrs2cam p2cam];
xdata = [Kvec nCamPts nLrsPts planePts lrsPts];
ydata = [camPts zeros(1,sum(nLrsPts))];
