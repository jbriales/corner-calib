classdef CCheckerboardComp < handle & CBaseComp
        
    properties
        Vasconcelos
        Zhang
    end
    
    methods
        
        % Constructor
        function obj = CCheckerboardComp ( )
            obj.Vasconcelos = CBaseComp;
            obj.Zhang = CBaseComp;
        end
        
        function optim( obj, checkerOptim )
            [R_v,t_v] = checkerOptim.optimizeRt_Vasc;
            obj.Vasconcelos.storeResult( [R_v t_v] );
            
            [R_z,t_z] = checkerOptim.optimizeRt_Zhang;            
            obj.Zhang.storeResult( [R_z t_z] );
        end
    end
           
end