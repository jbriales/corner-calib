% USES
%   [Tlrs2cam, Tp2cam, Err, Err0] = optimizeLaserCamCalib(Tlrs2cam,Tp2cam,Points,CamCalib)
%   
%   [Tlrs2cam, Err, Err0] = optimizeLaserCamCalib(Tlrs2cam,Tp2cam,Points,CamCalib)

function [Tlrs2cam, Tp2cam, Err, Err0] = optimizeLaserCamCalib(Tlrs2cam,Tp2cam,Points,CamCalib)

% refine calibration and plane poses    
if nargout == 4
    [x0,xdata,ydata] = laserCamOptData(Points,CamCalib,Tp2cam,Tlrs2cam);
    y                = reprojectLaserCamCalib(x0,xdata);
    Err0  = y-ydata;
    
    % levenberg-marquadt error minimization
    x = lsqcurvefit(@reprojectLaserCamCalib,x0,xdata,ydata);
    
    nPlanes = size(Tp2cam,3);
    
    for i=1:nPlanes
        
        rp2cam = [x(6*i+1) x(6*i+2) x(6*i+3)];
        tp2cam = [x(6*i+4) x(6*i+5) x(6*i+6)];
        Tp2cam(:,:,i) = [aa2R(rp2cam) tp2cam.'; 0 0 0 1];
    end
    
    rlrs2cam = [x(1); x(2); x(3)];
    tlrs2cam = [x(4); x(5); x(6)];
    Tlrs2cam = [aa2R(rlrs2cam) tlrs2cam; 0 0 0 1];
    
    % Reprojection error
    y   = reprojectLaserCamCalib(x,xdata);
    Err = y-ydata;

% refine calibration only
else
    [x0,xdata,ydata] = laserOptData(Points,Tp2cam,Tlrs2cam);    
    y                = reprojectLaserCalib(x0,xdata);
    Err0  = y-ydata;
    
    % levenberg-marquadt error minimization
    x = lsqcurvefit(@reprojectLaserCalib,x0,xdata,ydata);
    
    rlrs2cam = [x(1); x(2); x(3)];
    tlrs2cam = [x(4); x(5); x(6)];
    Tlrs2cam = [aa2R(rlrs2cam) tlrs2cam; 0 0 0 1];
    
    % Reprojection error
    y   = reprojectLaserCalib(x,xdata);
    Err = y-ydata;
    
    Tp2cam = Err;
    Err   = Err0;
end

    