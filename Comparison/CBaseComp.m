classdef CBaseComp < handle
    
    % Define static variables through static methods for all common indexes
    methods (Static = true)
        function val = idx_cam_sd(newval)
            persistent currentval;
            if nargin >= 1
                currentval = newval;
            end
            val = currentval;
        end
        
        function val = idx_scan_sd(newval)
            persistent currentval;
            if nargin >= 1
                currentval = newval;
            end
            val = currentval;
        end
        
        function val = idx_N_co(newval)
            persistent currentval;
            if nargin >= 1
                currentval = newval;
            end
            val = currentval;
        end
        
        function val = idx_Nsim(newval)
            persistent currentval;
            if nargin >= 1
                currentval = newval;
            end
            val = currentval;
        end
    end
    
end