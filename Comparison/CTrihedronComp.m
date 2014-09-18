classdef CTrihedronComp < handle & CBaseComp
    
    properties
        %Linear
        %DiagWeighted
        NonWeighted3D        
        Weighted3D
        NonWeighted2D
        Weighted2D
        Global
    end
    
    methods
        
        % Constructor
        function obj = CTrihedronComp ( )
            %obj.Linear          = CBaseComp;
            %obj.DiagWeighted    = CBaseComp;
            obj.NonWeighted3D    = CBaseComp;
            obj.Weighted3D       = CBaseComp;
            obj.NonWeighted2D    = CBaseComp;
            obj.Weighted2D       = CBaseComp;
            obj.Global           = CBaseComp;
        end
        
        function optim( obj, triOptim, WITHRANSAC )
            if WITHRANSAC
                triOptim.filterRotationRANSAC;
                triOptim.disp_N_R_inliers;
            end            
            %R_c_s_dw = triOptim.optimizeRotation_DiagWeighted;
            R_c_s_nw = triOptim.optimizeRotation_NonWeighted;            
            R_c_s_w  = triOptim.optimizeRotation_Weighted;
            
            if WITHRANSAC
                triOptim.filterTranslationRANSAC( R_c_s_w ); % Should receive some estimated rotation
                triOptim.disp_N_t_inliers;
            end
            R0_for_t = R_c_s_w;
            
            t_3D_nw = triOptim.optimizeTranslation_3D_NonWeighted( R0_for_t );
            t_3D_w  = triOptim.optimizeTranslation_3D_Weighted( R0_for_t );
            
            t_2D_nw = triOptim.optimizeTranslation_2D_NonWeighted( R0_for_t );
            t_2D_w = triOptim.optimizeTranslation_2D_Weighted( R0_for_t );

            [R_global, t_global] = triOptim.optimizeGlobal_Ort_3D( R_c_s_w, t_3D_w );  
            
            % Indexes to store are static variables in class CBaseComp
            %obj.DiagWeighted(cam_sd_N_it, scan_sd_N_it, Nsim_it) = {triOptim.optimizeRotation_DiagWeighted};
            obj.NonWeighted3D.storeResult( [ R_c_s_nw t_3D_nw] );
            obj.Weighted3D.storeResult( [ R_c_s_w  t_3D_w] );          
            obj.NonWeighted2D.storeResult( [ R_c_s_nw t_2D_nw] );
            obj.Weighted2D.storeResult( [ R_c_s_w  t_2D_w] );
            obj.Global.storeResult( [ R_global t_global] );
            
            if WITHRANSAC
                triOptim.resetRotationRANSAC();
                triOptim.resetTranslationRANSAC();
            end
            
        end
        
    end
            
end
            