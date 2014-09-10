classdef CCornerObservation
    %CCornerObservation Class to store data extracted from a corner frame
    %grabbed from a Rig
    %   Detailed explanation goes here
    
    properties
        % Information seen from Camera (for rotation) 
        R_c_w       % 3(3x1) Normals of World planes X,Y,Z
        A_R_c_w     % Uncertainty of set of 3 normals (9x9 matrix, rank 3)
        A_eps_c_w   % Minimal representation of uncertainty (3x3 matrix)
        
        % Information seen from Lidar (for rotation)
        v   % 3(2x1) 2D Direction of Lidar segment
        A_v % Minimal representation of uncertainty of set of 3 directions (3x3 matrix)
        A_lh    % (*) Uncertainty of 2D Homogeneous lines in Lidar
        
        % Information seen from Camera (for translation)
        l   % 3(3x1) 2D homogeneous lines in image
        
        % Information seen from Lidar (for translation)
        q   % 3(2x1) 2D intersection points of Lidar segments
        A_q % Uncertainty in 2D intersection points (2x2 matrix)
    end
    
    methods
    end
    
    %% Store data for final optimisation
            %tic
            lab = [1 2 3];
            co.complete = CompleteCO;
            co.thereis_line = thereis_line;        % Line - rotation correspondence
            co.lab_line = lab(thereis_line);
            co.thereis_corner = thereis_corner;    % Corner - translation correspondence
            co.lab_corner = lab(thereis_corner);
            

            
            
            % Do whatever necessary here
%             corresp = cell(1,3);
            corresp = co;
    
end

