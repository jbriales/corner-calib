% x = [rlrs2cam tlrs2cam rp2cam1 tp2cam1 ... rp2camN tp2camN]
% xdata = [Kvec nPlanePts nLrsPts planePts lrsPts]
% y = [camPts1 ... camPtsN 0 0 ... 0]

function y = reprojectLaserCalib(x,xdata)

rlrs2cam = [x(1); x(2); x(3)];
tlrs2cam = [x(4); x(5); x(6)];
Rlrs2cam = aa2R(rlrs2cam);

nPlanes = xdata(1);
nLrsPts = xdata(2:nPlanes+1);
p2cam   = xdata((nPlanes+2):(1+7*nPlanes));


lrsPtsIdx   = 2+7*nPlanes;

% num_data = (length(xdata)-lrsPtsIdx+1)/3;
% ylrs = [];
% tic
err = cell(1,nPlanes);
for i=1:nPlanes  
    clear L
    
    rp2cam = p2cam((6*i-5):(6*i-3)).';
    tp2cam = p2cam((6*i-2):(6*i)).';  
    Rp2cam = aa2R(rp2cam);
    
    % Get laser points
    LrsPts    = reshape(xdata(lrsPtsIdx:lrsPtsIdx+3*nLrsPts(i)-1),3,[]);
    lrsPtsIdx = lrsPtsIdx + 3*nLrsPts(i);      
    
    % Get Plane in Laser Frame
    PIlrs = [Rlrs2cam [0;0;0]; -tlrs2cam.'*Rlrs2cam 1]*[Rp2cam [0;0;0]; -tp2cam.'*Rp2cam 1]*[0; 0; 1; 0];
    % Get Laser Lines
    L(1:3,:) = [LrsPts(1:2,:); zeros(1,size(LrsPts,2))];
%     for j=1:size(LrsPts,2)
%         L(4:6,j) = cross(LrsPts(:,j),[0;0;10]);
%     end
    % Improve speed
    L(4:6,:) = -skew([0;0;10]) * LrsPts;
    % Get intersection between lines and plane
%     P=IntersectionLinePlane(PIlrs, L);
    % Avoid call function
    P = makeinhomogeneous( [PIlrs(4)*eye(3), skew_symetric_v(PIlrs(1:3));
                               -PIlrs(1:3)',                   zeros(1,3)] * L );
    % Get error between intersections and laser points
    err{i} = sqrt(sum((LrsPts(1:2,:,1)-P(1:2,:)).^2));
%     err    = sqrt(sum((LrsPts(1:2,:,1)-P(1:2,:)).^2));
%     ylrs   = [ylrs err];
end
ylrs = cell2mat( err );
% toc

y = ylrs;

end

function X = skew(x)
% Creates a skew symmetrix matrix from a vectr x in R3

X = [0 -x(3) x(2) ; x(3) 0 -x(1) ; -x(2) x(1) 0 ];
end
