% s7_checkResults

% Load again a set of lines
L = load(fullfile(path,'lines',file.lines)); 
h = figure();

imshow(Ic), hold on

% Plot detected lines
vx = 1:2*res_x;
vy = 1:2*res_y;
for i=1:3
    l = L(i,1:3);
    l(3) = l(3) - l(1) - l(2); % Correct 1-pixel of OpenCV to Matlab indexes
    plotHomLine(vx,vy,l,'r');
end

% Plot optimized intersection of lines
im_o = makeinhomogeneous(K * o);
plot( im_o(1), im_o(2), 'r*' )
clear im_o % Temporal variable for o position in pixels

% Plot current projection of scan onto image (with GT):
% ---------
% TODO: include translation
Npts = size(pts,2);
pts_proj = R_lev(:,1:2) * pts + repmat(t_lin,[1 Npts]);
pts_proj = K * pts_proj;
pts_proj = makeinhomogeneous( pts_proj );
pts_proj( :, pts_proj(1,:)<=0 | pts_proj(2,:)<=0 | pts_proj(1,:)>res_x | pts_proj(2,:)>res_y ) = [];
figure(h), hold on
plot(pts_proj(1,:),pts_proj(2,:),'g.')
% ---------
% END Plot

% Plot current projection of scan onto image (with linear estimation of
% translation):
% ---------
% TODO: include translation
Npts = size(pts,2);
pts_proj = R_c_s(:,1:2) * pts + repmat(t_c_s,[1 Npts]);
pts_proj = K * pts_proj;
pts_proj = makeinhomogeneous( pts_proj );
pts_proj( :, pts_proj(1,:)<=0 | pts_proj(2,:)<=0 | pts_proj(1,:)>res_x | pts_proj(2,:)>res_y ) = [];
figure(h), hold on
plot(pts_proj(1,:),pts_proj(2,:),'go')
% ---------
% END Plot

% Plot projection of LRF corners with GT Rt
lab = 'xyz';
for i=1:3
    ax = lab(i);

    c_proj = R_c_s(:,1:2) * c.(ax) + t_c_s;
    c_proj = K * c_proj;
    c_proj = makeinhomogeneous( c_proj );
    plot(c_proj(1),c_proj(2),'y*','LineWidth',5)
%     fprintf('GT projection axis %c:\n',ax)
%     disp(c_proj)
end

% Plot projection of LRF corners with linear estimated Rt
lab = 'xyz';
for i=1:3
    ax = lab(i);

    c_proj = R_lev(:,1:2) * c.(ax) + t_lin;
    c_proj = K * c_proj;
    c_proj = makeinhomogeneous( c_proj );
    plot(c_proj(1),c_proj(2),'c*','LineWidth',5)
%     fprintf('Estimated projection axis %c:\n',ax)
%     disp(c_proj)
end

% Plot projection of LRF corners with linear estimated Rt
lab = 'xyz';
for i=1:3
    ax = lab(i);

    c_proj = R_lev(:,1:2) * c.(ax) + t_img;
    c_proj = K * c_proj;
    c_proj = makeinhomogeneous( c_proj );
    plot(c_proj(1),c_proj(2),'m*','LineWidth',5)
%     fprintf('Estimated projection axis %c:\n',ax)
%     disp(c_proj)
end

hold off