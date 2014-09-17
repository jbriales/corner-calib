classdef CComparison < handle
    %CComparison class to grab data from several simulations(TODO: comment)
    
    % Constructor is object-dependent

    properties 
        TrihedronComp   % Trihedron method
        KwakComp        % Kwak's method
        ZhangComp       % Zhang's algorithm
        VasconcelosComp % Vasconcelos' method
        
        N_sim           % Number of simulations performed for each tuple
        
        cam_sd_n        % Camera noise levels
        scan_sd_n       % Lidar noise levels
        N_co_n          % Number of correspondences levels
        
        Rt0              % Store groundtruth Rt for this specific comparison
    end
    
    methods
               
        % Constructor TODO: add propertie GTComp if the GT change every it
        function obj = CComparison( Rt0, N_sim, cam_sd_n, scan_sd_n, N_co_n )
            
            cam_sd_N  = size(cam_sd_n,2);
            scan_sd_N = size(scan_sd_n,2);
            N_co_N    = size(N_co_n,2);            
            
            obj.TrihedronComp   = CTrihedronComp(N_sim, cam_sd_N, scan_sd_N, N_co_N);     
            obj.KwakComp        = CKwakComp(N_sim, cam_sd_N, scan_sd_N, N_co_N);          
            obj.ZhangComp       = CZhangComp(N_sim, cam_sd_N, scan_sd_N, N_co_N);         
            obj.VasconcelosComp = CVasconcelosComp(N_sim, cam_sd_N, scan_sd_N, N_co_N);   
            
            obj.Rt0 = Rt0; % Store groundtruth transformation
            
            obj.cam_sd_n        = cam_sd_n;
            obj.scan_sd_n       = scan_sd_n;
            obj.N_co_n          = N_co_n;
            
            obj.N_sim           = N_sim;
        end
        
        function dispAvailable( obj )
            fprintf('==========================================\n');
            fprintf('Object data:\n')
            fprintf('==========================================\n');
            fprintf('Number of simulations for each tuple: %d\n', obj.N_sim);
            fprintf('Vector of camera noises (pixels):\n');
            disp(obj.cam_sd_n);
            fprintf('Vector of LRF noises (m):\n');
            disp(obj.scan_sd_n);
            fprintf('Vector of N of observations:\n');
            disp(obj.N_co_n);
            fprintf('==========================================\n\n');
        end
        
        
               
        function obj = plotCameraNoise( obj )
            % Set GT from object properties
            x_gt = obj.Rt0;
            
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

            % Calculate the errors
            R_gt = x_gt(1:3,1:3); t_gt = x_gt(1:3,4);
            cam_sd_vec = '';
            for i = 1:size(obj.cam_sd_n,2)
                cam_sd = obj.cam_sd_n(i);
                cam_sd_vec = [cam_sd_vec {num2str(cam_sd)}];
                for j = 1:obj.N_sim

                    % Trihedron error
                    if WITHTRIHEDRON
                        x = obj.TrihedronComp.(TRIHEDRON_METHOD){i,scan_sd_it,N_co_it,j};
                        R_err_trih(i,j) = angularDistance(x(1:3,1:3),R_gt);
                        t_err_trih(i,j) = norm(x(1:3,4)-t_gt);
                    end
                    
                    if WITHWASIEL
                        x_k = obj.KwakComp.Wasielewski{i,scan_sd_it,N_co_it,j};
                        R_err_wasiel(i,j) = angularDistance(x_k(1:3,1:3),R_gt);
                        t_err_wasiel(i,j) = norm(x_k(1:3,4)-t_gt);
                    end
                    
                    % Corner error
                    if WITHKWAK
                        x_k = obj.KwakComp.(KWAK_METHOD){i,scan_sd_it,N_co_it,j};
                        R_err_kwak(i,j) = angularDistance(x_k(1:3,1:3),R_gt);
                        t_err_kwak(i,j) = norm(x_k(1:3,4)-t_gt);
                    end

                    % Zhang error
                    if WITHZHANG
                        x_z = obj.ZhangComp.Linear{i,scan_sd_it,N_co_it,j};
                        R_err_zhang(i,j) = angularDistance(x_z(1:3,1:3),R_gt);
                        t_err_zhang(i,j) = norm(x_z(1:3,4)-t_gt);
                    end

                    % Vasconcelos error
                    if WITHVASC
                        x_v = obj.VasconcelosComp.Linear{cam_sd,scan_sd_it,N_co_it,j};
                        R_err_vasc(i,j) = angularDistance(x_v(1:3,1:3),R_gt);
                        t_err_vasc(i,j) = norm(x_v(1:3,4)-t_gt);
                    end
                end
            end
            
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
        
        function obj = plotLidarNoise( obj, x_GT )
            %TODO            
        end
        
    end
    
end