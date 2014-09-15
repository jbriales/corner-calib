classdef CCornerObservation
    %CCornerObservation Class for storage of Camera and LRF data in a
    %frame
    %   Detailed explanation goes here
    
    properties
        % Camera properties
        cam_l       % 3(3x1) cell array of homogeneous lines in image (coincident with direction of normals to reprojection planes from camera center after normalization)
        
        % LRF properties
        LRF_q       % 3(2x1) cell array of 2D intersection points of scan lines
    end
    
    methods
        function obj = CCornerObservation( cam_l, LRF_q )
            obj.cam_l = cam_l;
            obj.LRF_q = LRF_q;
        end
        
    end
    
end
