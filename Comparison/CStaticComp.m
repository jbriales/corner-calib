classdef CStaticComp < handle
    % CStaticComp Parent class with static variables and methods for
    % controlling indexing of all simulations with common parameters
    
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
        
        function val = idx_N_sim(newval)
            persistent currentval;
            if nargin >= 1
                currentval = newval;
            end
            val = currentval;
        end
    end
    
    % Define static variables through static methods for all dimensions
    methods (Static = true)
        function val = dim_cam_sd(newval)
            persistent currentval;
            if nargin >= 1
                currentval = newval;
            end
            val = currentval;
        end
        
        function val = dim_scan_sd(newval)
            persistent currentval;
            if nargin >= 1
                currentval = newval;
            end
            val = currentval;
        end
        
        function val = dim_N_co(newval)
            persistent currentval;
            if nargin >= 1
                currentval = newval;
            end
            val = currentval;
        end
        
        function val = dim_N_sim(newval)
            persistent currentval;
            if nargin >= 1
                currentval = newval;
            end
            val = currentval;
        end
    end
    
end