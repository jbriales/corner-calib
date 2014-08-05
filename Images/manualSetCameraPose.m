function [R_w_c, t_w_c] = manualSetCameraPose( t_w_c, angle )

hF = figure; hold on
axis equal
rotate3d on
axis([-1 +1 -1 +1 -1 +1])
plotframe(eye(4),1,'W','k')

% User input
if ~exist('t_w_c','var')
    t_w_c = [];
    while ~( size(t_w_c,1)==3 && size(t_w_c,2)==1 )
        str = input('Camera position wrt world: ','s');
        if ~isempty(str)
            t_w_c = textscan( str, '%f' );
            t_w_c = t_w_c{:};
        else
            t_w_c = [1 1 -1]';
        end
    end
    
    [az,el] = view;
end
[az,el] = view;

% Set approximate rotation (with X in XY plane)
Rz = - t_w_c / norm(t_w_c);
Rx = cross( Rz, [0 0 1]');
Rx = Rx / norm(Rx);
Ry = cross( Rz, Rx );
R_w_c = [Rx Ry Rz];

plotframe([ R_w_c t_w_c ; 0 0 0 1 ], 0.5, 'c', 'b');

if exist('angle','var')
    R_w_c = R_w_c * RotationZ(deg2rad(angle));
else
    disp('Set camera rotation')
    disp(' exit - keeps current rotation and exit')
    disp(' rot  - change current rotation')
    while 1
        control = input('Choose action: ','s');
        switch control
            case 'rot'
                angle = input('Angle: ');
                Rot = RotationZ(deg2rad(angle));
                
                clf
                subplot(121)
                hold on
                axis equal
                rotate3d on
                axis([-1 +1 -1 +1 -1 +1])
                plotframe(eye(4),1,'W','k')
                view(az,el)
                plotframe([ R_w_c * Rot t_w_c ; 0 0 0 1 ], 0.5, 'c', 'b');
                
                subplot(122)
                plotEstimatedImg( (R_w_c * Rot)' ) % Input is R_c_w so it is transposed
            case 'exit'
                break
        end
    end
    R_w_c = R_w_c * Rot;
end

close(hF);

end

function plotEstimatedImg( R_c_w )

title('Approximate image with current estimated rotation in Camera')
hold on
rgb = 'rgb';
for i=1:3
    p  = makeinhomogeneous(R_c_w(:,i));
    op = [ zeros(2,1) p ];
    plot( op(1,:), op(2,:), rgb(i) )
end
axis equal

end