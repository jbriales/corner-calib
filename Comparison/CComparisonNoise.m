classdef CComparisonNoise < handle
    %CComparison class to grab data from several simulations(TODO: comment)
    
    % Constructor is object-dependent

    properties 
        TrihedronComp   % Trihedron method
        KwakComp        % Kwak's method
        ZhangComp       % Zhang's algorithm
        VasconcelosComp % Vasconcelos' method
        N_samples
    end
    
    methods
               
        % Constructor 
        function obj = CComparisonNoise( N_sim, cam_sd_N, scan_sd_N, N_samples )
        
            obj.TrihedronComp   = CTrihedronComp(N_sim, cam_sd_N, scan_sd_N);     %(TODO)
            obj.KwakComp        = CKwakComp(N_sim, cam_sd_N, scan_sd_N);          %(TODO)
            obj.ZhangComp       = CZhangComp(N_sim, cam_sd_N, scan_sd_N);         %(TODO)
            obj.VasconcelosComp = CVasconcelosComp(N_sim, cam_sd_N, scan_sd_N);   %(TODO)
            obj.N_samples       = N_samples;
            
        end
        
    end
    
end