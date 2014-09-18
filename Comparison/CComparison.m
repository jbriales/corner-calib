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
            fprintf('==========================================\n\n');
        end
        function patterns = getPatterns( obj )
            patterns = fieldnames( obj );
        end
        function methods = getMethods( obj, pattern )
            methods = fieldnames( obj.(pattern) );
            methods(end) = []; % Remove mem
        end
        
        function obj = plotCameraNoise( obj )
            % Set GT from object properties
            x_gt = obj.Rt_gt;
            
            % Plot options
            plot_sim_file = fullfile( pwd, 'plotCameraNoise.ini' );
            plotOpts = readConfigFile( plot_sim_file );
            extractStructFields( plotOpts );
            clear plotOpts;

            % Check if the scan_sd value is correct
            if scan_sd_it > size(obj.scan_sd_n,2)
                error('Value of "scan_sd_it" out of range (max = %i)',size(obj.scan_sd_n,2));
            end
            scan_sd = obj.scan_sd_n(scan_sd_it);    

            % Check if the N_samples value is correct 
            if N_co_it > size(obj.N_co_n,2)
                error('Value of "N_co_it" out of range (max = %i)',size(obj.N_co_n,2));
            end
            N_co    = obj.N_co_n(N_co_it);

            % Extract dim x Nsim matrices for representation
            all_R_err = [];
            all_t_err = [];
            for idx_pattern = 1:size(patterns,1)
                pat = patterns{idx_pattern,1};
                methods = patterns{idx_pattern,2};
                for idx_method = 1:length(methods)
                    met = methods{idx_method};
                    Rt = squeeze(obj.(pat).(met).mem(:,scan_sd_it,N_co_it,:));
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
                end
            end
            err.xtick = num2str( obj.cam_sd_n );
            
            % TODO: Complete options (color, grouping, etc.)
            boxplot(R_err,'plotstyle','compact');
            
            % Plot the errors
            if WITHTRIHEDRON
                R_err = [R_err_trih'];
                t_err = [t_err_trih'];
                b = [repmat({'Trihedron'},1,size(obj.cam_sd_n,2))];
            end
            if WITHWASIEL
                R_err = [R_err R_err_wasiel'];
                t_err = [t_err t_err_wasiel'];
                b = [b, repmat({'Wasiel'},1,size(obj.cam_sd_n,2))];
            end
            if WITHKWAK
                R_err = [R_err R_err_kwak'];
                t_err = [t_err t_err_kwak'];
                b = [b, repmat({'Kwak'},1,size(obj.cam_sd_n,2))];
            end
            if WITHZHANG
                R_err = [R_err R_err_zhang'];
                t_err = [t_err t_err_zhang']; 
                b = [b, repmat({'Zhang'},1,size(obj.cam_sd_n,2))];
            end
            if WITHVASC
                R_err = [R_err R_err_vasc'];
                t_err = [t_err t_err_vasc'];
                b = [b, repmat({'Vasconcelos'},1,size(obj.cam_sd_n,2))];
            end
            
            N_plots = WITHTRIHEDRON + WITHWASIEL + WITHKWAK + WITHZHANG + WITHVASC;
            color   = color(1:N_plots,:);                           
            a = repmat(cam_sd_vec,1,N_plots);
%             b = [repmat({'Trihedron'},1,3), repmat({'Kwak'},1,size(obj.cam_sd_n,2))];
            figure
            subplot(211)
            boxplot(R_err,{a,b},'colors', repmat(color,size(obj.cam_sd_n,2),1), 'factorgap',[5 0.05],'plotstyle','compact');
            set(gca,'YScale','log') 
            subplot(212)
            boxplot(t_err,{a,b},'colors', repmat(color,size(obj.cam_sd_n,2),1), 'factorgap',[5 0.05],'plotstyle','compact');
            set(gca,'YScale','log')
%             set(gca,'xticklabel',{'Direct care','Housekeeping','Mealtimes','Medication','Miscellaneous','Personal care'})
        end
        
        function [R_err] = R_err( obj, Rt )
            R_err = angularDistance( Rt(1:3,1:3), obj.Rt_gt(1:3,1:3) );
        end
        function [t_err] = t_err( obj, Rt )
            t_err = norm( Rt(1:3,4) - obj.Rt_gt(1:3,4) );
        end
        
        function obj = plotLidarNoise( obj, x_GT )
            %TODO            
        end
        
    end
    
end