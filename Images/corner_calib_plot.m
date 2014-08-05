function corner_calib_plot(x, x0, img_gray, pts, w)

    % Function that plots in the complete image the solution, the initial
    % estimation and the rectangle of the subimage
    
    % Inputs:
    % x         - final estimation [px py t1 t2 t3]'
    % x0        - initial estimation [px py t1 t2 t3]'
    % img_gray  - complete gray image

    hold off
    imshow(img_gray)
    hold on
    
    mh = @(x) [ x ; 1 ];
        
    % Old line
    if ~isempty(x0)
        x0 = x0(:); % Assure column vector
        p = x0(1:2);
        ang = x0(3:5);
        v = [ cos( ang ), sin( ang ) ]';
        rgb = 'rgb';
        for k=1:3
            hom_line = cross(mh(p),mh(p+v(:,k)));
            plotHomLineWin( hom_line, ['-.',rgb(k)] );
        end
    end
    
    % New line
    x = x(:); % Assure column vector
    p = x(1:2);
    ang = x(3:5);
    v = [ cos( ang ), sin( ang ) ]';
    rgb = 'rgb';
    for k=1:3
        hom_line = cross(mh(p),mh(p+v(:,k)));
        plotHomLineWin( hom_line, rgb(k) );
    end
    
    if exist('pts','var') && exist('w','var')
        hFun1 = plotImgFun( pts{1}, w{1} );
        hFun2 = plotImgFun( pts{2}, w{2} );
        hFun3 = plotImgFun( pts{3}, w{3} );
        set(hFun1,'Visible','off')
        set(hFun2,'Visible','off')
        set(hFun3,'Visible','off')
    end
end

