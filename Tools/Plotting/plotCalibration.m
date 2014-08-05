function plotCalibration( imgs, scans, gt )

N = length( imgs );

hF = figure; hold on
axis equal
rotate3d on

title('Calibration scene')

% World frame
plotframe( eye(4), 5, 'W', 'k' )

t = [scans.t_w_s];
x_s = t(1,:);
y_s = t(2,:);
z_s = t(3,:);
for i=1:N
    hTraj = plot3(x_s,y_s,z_s,'-b');
    
    R_w_s = scans(i).R_w_s;
    t_w_s = scans(i).t_w_s;
    hFr = plotframe( [R_w_s t_w_s ; 0 0 0 1], 0.5, 'S', 'g' );
    
    pts = scans(i).xy;
    pts = R_w_s(:,1:2) * pts + repmat(t_w_s,1,length(pts));
    hPts = plot3( pts(1,:), pts(2,:), pts(3,:), '.b' );
    
    figure(hF)
    pause
    delete(hTraj)
    delete(hFr)
    delete(hPts)
end