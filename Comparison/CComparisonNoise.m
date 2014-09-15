classdef CComparisonNoise < handle
    %CComparison class to grab data from several simulations(TODO: comment)
    
    % Constructor is object-dependent

    properties 
        TrihedronComp   % Trihedron method
        KwakComp        % Kwak's method
        ZhangComp       % Zhang's algorithm
        VasconcelosComp % Vasconcelos' method
        N_samples       % Number of samples per calibration
        N_sim           % Number of simulations performed
        cam_sd_N        % Number of camera noise levels
        scan_sd_N       % Number of lidar noise levels
    end
    
    methods
               
        % Constructor TODO: add propertie GTComp if the GT change every it
        function obj = CComparisonNoise( N_sim, cam_sd_N, scan_sd_N, N_samples )
        
            obj.TrihedronComp   = CTrihedronComp(N_sim, cam_sd_N, scan_sd_N);     
            obj.KwakComp        = CKwakComp(N_sim, cam_sd_N, scan_sd_N);          %(TODO)
            obj.ZhangComp       = CZhangComp(N_sim, cam_sd_N, scan_sd_N);         %(TODO)
            obj.VasconcelosComp = CVasconcelosComp(N_sim, cam_sd_N, scan_sd_N);   %(TODO)
            obj.N_samples       = N_samples;
            
        end
        
        function obj = plotCameraNoise( obj, WITHTRIHEDRON, WITHCORNER, WITHZHANG )
            %TODO            
        end
        
        function obj = plotLidarNoise( obj, WITHTRIHEDRON, WITHCORNER, WITHZHANG )
            %TODO            
        end
        
    end
    
end