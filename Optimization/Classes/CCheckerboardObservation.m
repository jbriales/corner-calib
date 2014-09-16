classdef CCheckerboardObservation
    %CCheckerboardObservation Class for storage of Camera and LRF data in a
    %frame
    %   Detailed explanation goes here
    
    properties
        % Grid properties
        plane_T      % 4x4 transformation matrix with pose of grid plane
        
        % LRF properties
        LRF_pts     % 2xN array of LRF points in [mm]!!!
    end
    
    methods
        function obj = CCheckerboardObservation( T, pts)
            obj.plane_T = T;
            obj.LRF_pts = pts;
        end
    end
    
end