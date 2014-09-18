classdef CComparison < handle & CStaticComp
    %CComparison class to grab data from several simulations(TODO: comment)
    
    % Constructor is object-dependent

    properties 
        Trihedron       % Trihedron methods
        Corner          % Corner methods
        Checkerboard    % Checkerboard methods
    end
    
    properties (Access = protected)
        cam_sd_n        % Camera noise levels
        scan_sd_n       % Lidar noise levels
        N_co_n          % Number of correspondences levels
        
        Rt_gt           % Store groundtruth Rt for this specific comparison
        N_sim
    end
    
    methods
               
        % Constructor TODO: add property GTComp if the GT change every it
        function obj = CComparison( Rt_gt, N_sim, cam_sd_n, scan_sd_n, N_co_n )
            
            obj.dim_cam_sd( size(cam_sd_n,2) );
            obj.dim_scan_sd( size(scan_sd_n,2) );
            obj.dim_N_co( size(N_co_n,2) );
            obj.dim_N_sim( N_sim );
            obj.N_sim = N_sim;
            
            obj.Trihedron    = CTrihedronComp;
            obj.Corner       = CCornerComp;
            obj.Checkerboard = CCheckerboardComp;
            
            obj.Rt_gt = Rt_gt; % Store groundtruth transformation
            
            obj.cam_sd_n        = cam_sd_n;
            obj.scan_sd_n       = scan_sd_n;
            obj.N_co_n          = N_co_n;
        end
        
        % Show object information
        function dispAvailable( obj )
            fprintf('==========================================\n');
            fprintf('Object data:\n')
            fprintf('==========================================\n');
            fprintf('Simulated methods:')
            patterns = obj.getPatterns;
            for i = 1:length(patterns)
                disp( patterns{i} );
                disp('-----------------');
                disp( obj.getMethods(patterns{i}) );
            end
            fprintf('Number of simulations for each tuple: %d\n', obj.N_sim);
            fprintf('Vector of camera noises (pixels):\n');
            disp(obj.cam_sd_n);
            fprintf('Vector of LRF noises (m):\n');
            disp(obj.scan_sd_n);
            fprintf('Vector of N of observations:\n');
            disp(obj.N_co_n);
            fprintf('Complete lines for .ini\n');
            fprintf('cam_sd_vals  = [ ');
            fprintf('%f ',obj.cam_sd_n); fprintf(']\n');
            fprintf('scan_sd_vals = [ ');
            fprintf('%f ',obj.scan_sd_n); fprintf(']\n');
            fprintf('Nobs_vals    = [ ');
            fprintf('%f ',obj.N_co_n); fprintf(']\n');
            fprintf('==========================================\n\n');
        end
        function patterns = getPatterns( obj )
            patterns = fieldnames( obj );
        end
        function methods = getMethods( obj, pattern )
            methods = fieldnames( obj.(pattern) );
            methods(end) = []; % Remove mem
        end
        
        function checkValue( obj, field, v )
            for val = v
                if isempty(find(val==obj.(field),1))
                    fprintf('%s available values:\n',field);
                    disp(obj.(field));
                    error('Non existent value %f for field %s',val,field);
                end
            end            
        end
        
        function indexes = getIndexes( obj, field, v )
            indexes = [];
            for val = v
                idx = find(val==obj.(field),1);
                if isempty(idx)
                    fprintf('%s available values:\n',field);
                    disp(obj.(field));
                    error('Non existent value %f for field %s',val,field);
                else
                    indexes(end+1) = idx;
                end
            end            
        end
        
        function obj = plotData( obj, cam_sd_vals, scan_sd_vals, Nobs_vals )           
            % Plot options
            plot_sim_file = fullfile( pwd, 'plotBoxplot.ini' );
            plotOpts = readConfigFile( plot_sim_file );
            extractStructFields( plotOpts );
            clear plotOpts;
            
            % Check size of inputs (only one should be > 1)
            s = [ length(cam_sd_vals), length(scan_sd_vals), length(Nobs_vals) ];
            if length(find(s>1)) > 1
                error('Currently only one variable can be a vector')
            end
            [Nx,imax] = max(s); % The number of groups (one for each X value)
            % Get label
            vals = {cam_sd_vals, scan_sd_vals, Nobs_vals};
            xlabels = {'Cam noise','LRF noise','Nobs'};
            val_label = {};
            for val = vals{imax};
                val_label{end+1} = num2str(val);
            end
            xlab = xlabels{imax};
            
            % Check if the scan_sd value is correct
            cam_idxs  = obj.getIndexes('cam_sd_n',cam_sd_vals);
            scan_idxs = obj.getIndexes('scan_sd_n',scan_sd_vals);
            Nobs_idxs = obj.getIndexes('N_co_n',Nobs_vals);

            % Extract dim x Nsim matrices for representation
            all_R_err = [];
            all_t_err = [];
            Cval = cell(1,0);
            Ctag = cell(1,0);
            N_met = 0;
            for idx_pattern = 1:size(patterns,1)
                pat = patterns{idx_pattern,1};
                methods = patterns{idx_pattern,2};
                for idx_method = 1:length(methods)
                    met = methods{idx_method};
%                     Rt = squeeze(obj.(pat).(met).mem(:,scan_sd_it,N_co_it,:));
                    Rt = squeeze(obj.(pat).(met).mem(cam_idxs,scan_idxs,Nobs_idxs,:));
                    if size(Rt,2)==1
                        warning('Check squeeze changed dimension order if too many singletons');
                        keyboard
                        Rt = Rt';
                    end
                    R_err = cell2mat( cellfun(@(x)obj.R_err(x), Rt, 'UniformOutput',false) );
                    t_err = cell2mat( cellfun(@(x)obj.t_err(x), Rt, 'UniformOutput',false) );
                    err.(pat).(met).R = R_err;
                    err.(pat).(met).t = t_err;
                    all_R_err = [ all_R_err , R_err' ];
                    all_t_err = [ all_t_err , t_err' ];
                    
                    met_label = repmat({met},1,Nx);
                    Cval = {Cval{:},val_label{:}};
                    Ctag = {Ctag{:},met_label{:}};
                    
                    N_met = N_met+1;
                end
            end
            err.xtick = num2str( obj.cam_sd_n );
            
            % TODO: Complete options (color, grouping, etc.)
            color = repmat(rand(N_met,3),Nx,1);
            Rlab = 'Rotation error (deg)';
            tlab = 'Translation error (m)';
            figure
            subplot(211)
%             boxplot(all_R_err,{Cval,Ctag},'colors', color, 'factorgap',[5 0.05],'plotstyle','compact');
            boxplot(all_R_err,{Cval,Ctag},'colors', color, 'factorgap',10,'plotstyle','compact');
            set(gca,'YScale','log')
            xlabel(xlab)
            ylabel(Rlab)
            subplot(212)
            boxplot(all_t_err,{Cval,Ctag},'colors', color, 'factorgap',[5 0.05],'plotstyle','compact');
            set(gca,'YScale','log')
            xlabel(xlab)
            ylabel(tlab)
%             set(gcf,'units','normalized','position',[0 0 1 1]);
        end
        
        function [R_err] = R_err( obj, Rt )
            R_err = angularDistance( Rt(1:3,1:3), obj.Rt_gt(1:3,1:3) );
        end
        function [t_err] = t_err( obj, Rt )
            t_err = norm( Rt(1:3,4) - obj.Rt_gt(1:3,4) );
        end
        
        function obj = plotCamNoise( obj )
            %TODO            
        end
        function obj = plotLidarNoise( obj )
            %TODO            
        end
        
        
        
    end
    
end