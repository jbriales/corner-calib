classdef CConfigCamera
    %CConfigCamera Config class for Camera object
    % This base class stores configuration parameters
    %   Constructor:
    %   config = CConfigCamera( K, res, f, sd  )
    %
    % Adjustable parameters are those inputs in constructor:
    %   K - calibration matrix
    %   res - image resolution (height x width)
    %   sd - Standard Deviation (in [pixel]) of image pixels
    
    properties (SetAccess = protected ) % Only changeable through constructor
        K       % Intrinsic calibration matrix
        res     % Image resolution (height x width)
        f       % Focal length
        sd      % Standard Deviation in image pixels
        
        ax      % Min and max values for image axes coordinates
        border  % Homogeneous lines for camera borders
        
        FOVh    % FOV horizontal in rads (x direction?)
        FOVv    % FOV vertical in rads (y direction?)
        FOVhd   % FOV horizontal in degs (x direction?)
        FOVvd   % FOV vertical in degs (y direction?)
    end
    
    properties (SetAccess = protected, Dependent)
        % Empty
    end
       
    methods
        % Constructor
        function obj = CConfigCamera( K, res,f, sd  )
            obj.K  = K;
            obj.res = res;
            obj.f  = f;
            obj.sd = sd;
            
            obj.ax     = [ 1 res(2) 1 res(1) ]; 
            obj.border = [ 1 1 0 0
                           0 0 1 1
                           obj.ax ];
                       
            % TODO: Check formulae
            obj.FOVh = 2 * atan( K(1,3)/K(1,1) );
            obj.FOVv = 2 * atan( K(2,3)/K(2,2) );
            obj.FOVhd = rad2deg( obj.FOVh );
            obj.FOVvd = rad2deg( obj.FOVv );
        end
    end
end