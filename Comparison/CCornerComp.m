classdef CCornerComp < handle & CBaseComp
    
    properties
        Kwak_NonWeighted
        Kwak_ConstWeighted
        Wasielewski
    end
    
    methods
        
        % Constructor
        function obj = CCornerComp ( )
            obj.Kwak_NonWeighted     = CBaseComp;
            obj.Kwak_ConstWeighted   = CBaseComp;
            obj.Wasielewski          = CBaseComp;
        end
        
        function optim( obj, cornerOptim )
            [R_k_nw, t_k_nw] = cornerOptim.optimizeRt_NonWeighted;
            [R_k_cw, t_k_cw] = cornerOptim.optimizeRt_ConstWeighted;
            [R_kC_nw, t_kC_nw] = cornerOptim.optimizeRt_C_NonWeighted;
                        
            obj.Kwak_NonWeighted.storeResult( [R_k_nw  t_k_nw] );
            obj.Kwak_ConstWeighted.storeResult( [R_k_cw  t_k_cw] );
            obj.Wasielewski.storeResult( [R_kC_nw t_kC_nw] );
        end
        
    end
        
end