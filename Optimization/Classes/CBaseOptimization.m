classdef CBaseOptimization < handle % handle for compatibility with C...Optimization
    %COptimization Base class with common properties for optimization
    %processes
    %   Detailed explanation goes here
    
    properties
        debug_level % Verbose level when optimizing
        maxIters    % Max number of iterations in LM optimization
        minParamChange   % Minimum change in param (norm2) to stop
        minErrorChange  % Minimum change in cost function to stop
    end
    
    methods
        function obj = CBaseOptimization( debug_level, maxIters, minParamChange, minErrorChange )
            if ~exist('debug_level','var')
                debug_level = 2;
            end
            obj.debug_level = debug_level;
            
            if ~exist('maxIters','var')
                maxIters = 50;
            end
            obj.maxIters = maxIters;
            
            if ~exist('minParamChange','var')
                minParamChange = 1e-8;
            end
            obj.minParamChange = minParamChange;
            
            if ~exist('minErrorChange','var')
                minErrorChange = 1e-8;
            end
            obj.minErrorChange = minErrorChange;
        end
        
        function x = optimize( Fun, x0, space, weighted )
            [ x, err, errNorm, W ] = LM_Man_optim(Fun,x0,...
                'space',space,'weighted',weighted,...
                'debug', obj.debug_level,...
                'maxIters', obj.maxIters,...
                'minChange', obj.minParamChange,...
                'minErrorChange', obj.minErrorChange);
        end
    end
    
end

