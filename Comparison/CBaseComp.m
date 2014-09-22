classdef CBaseComp < handle & CStaticComp
    % CBaseComp Basic object for storage of results and some elementary
    % operations with that data common for all the methods
    
    properties
%     properties (Access = protected)
        mem
        cov_R
        cov_t
        cov_Rt
    end
    
    methods
        % Constructor
        function obj = CBaseComp( )
            obj.mem    = obj.preallocate_cell;
            obj.cov_R  = obj.preallocate_cell;
            obj.cov_t  = obj.preallocate_cell;
            obj.cov_Rt = obj.preallocate_cell;
        end
        
        % Method to store information in current indexes
        function obj = storeResult( obj, input )
            obj.mem{ obj.idx_cam_sd,...
                     obj.idx_scan_sd,...
                     obj.idx_N_co,...
                     obj.idx_N_sim } = input;
        end
        function obj = store( obj, field, input )
            obj.(field){ obj.idx_cam_sd,...
                         obj.idx_scan_sd,...
                         obj.idx_N_co,...
                         obj.idx_N_sim } = input;
        end
        
        % Method to extract information
%         function arr = extractDimension( obj, dim )
%             switch dim
%                 case 'cam_sd'
%                     arr = obj.mem( :, 1, 1, : );
%                 case 'scan_sd'
%                     arr = obj.mem( 1, :, 1, : );
%                 case 'N_co'
%                     arr = obj.mem( 1, 1, :, : );
%                 otherwise
%                     warning('Non valid dimension tag')
%                     arr = {};
%             end
%         end
        
        % Auxiliar
        function empty_cell = preallocate_cell( obj )
            empty_cell = cell(obj.dim_cam_sd,...
                              obj.dim_scan_sd,...
                              obj.dim_N_co,...
                              obj.dim_N_sim);
        end
    end
    
end