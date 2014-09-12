function [x0,xdata,ydata] = laserOptData(Points,Tp2cam,Tlrs2cam)

rlrs2cam = vrrotmat2vec(Tlrs2cam(1:3,1:3));
rlrs2cam = rlrs2cam(1:3)*rlrs2cam(4);
tlrs2cam = Tlrs2cam(1:3,4);
lrs2cam  = [rlrs2cam tlrs2cam.'];

p2cam    = [];
lrsPts   = [];
nPlanes  = size(Tp2cam,3); 

for i=1:nPlanes
    rp2cam = vrrotmat2vec(Tp2cam(1:3,1:3,i));
    rp2cam = rp2cam(1:3)*rp2cam(4);
    tp2cam = Tp2cam(1:3,4,i);
    p2cam  = [p2cam rp2cam tp2cam.'];    
    lrsPts     = [lrsPts reshape([Points{i}(1:2,:); zeros(1,size(Points{i},2))],1,[])];
    nLrsPts(i) = size(Points{i},2);
end


x0    = lrs2cam;
xdata = [nPlanes nLrsPts p2cam lrsPts];
ydata = zeros(1,sum(nLrsPts));
