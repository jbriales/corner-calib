classdef COdometryTest < handle
    %COdometryTest Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        WITH_ELQURSH = true;
        WITH_TRIHEDRON = true;
        WITH_PLOT = true;
        WITH_COVARIANCE = true;
        cam_sd = 0;
        lin_sd = 1e-5;
        cube_sd = 1e-3;
        rig_ini = 'rig.ini';
        pose_gen_ini = 'pose_gen.ini';
        
        % Stored variables
        R_gt
        
        % Temporary objects
        prim2D % 2D primitive
        
        % Graphic handles
        hFigure
        hF_scene
        hF_image
        
        % Used objects
        poseFactory
        Cam
        Cube
    end
    
    methods
        function this = COdometryTest( WITH_PLOT,...
                cam_sd, cube_sd,...
                rig_ini, pose_gen_ini )
            
            setConstructorInputs;
            
            % Read camera  options
            rig_config_file = fullfile( pwd, this.rig_ini );
            rigOpts = readConfigFile( rig_config_file );
            extractStructFields( rigOpts );
            
            % Generate random pose for camera
            gen_conf_F = fullfile( pwd, this.pose_gen_ini );
            this.poseFactory = CRandomPoses( readConfigFile( gen_conf_F,'[Vanishing]' ) );
            [R_w_c, t_w_c] = this.poseFactory.gen( 1 );
            this.R_gt = R_w_c';
            
            % Create simulated camera
            this.Cam = CSimCamera( R_w_c, t_w_c, K, res, f, this.cam_sd );
            
            % Create pattern
            this.Cube = CHalfCube( 1, eye(3), zeros(3,1), this.cube_sd );
            
            if this.WITH_PLOT
                this.hFigure  = figure;
                this.hF_scene = subplot(1,2,1);
                this.Cube.plotScene(this.Cam);
                
                this.hF_image = subplot(1,2,2); hold on;
            end
        end
        
        function set.WITH_PLOT( this, value )
            this.WITH_PLOT = value;
            if ( value ) % Activation case
                this.hFigure  = figure;
                this.hF_scene = subplot(1,2,1);
                this.Cube.plotScene(this.Cam);
                
                this.hF_image = subplot(1,2,2); hold on;
            else % Deactivation case
                close all
            end
        end
        function set.cam_sd( this, value )
            if value < 0
                warning('Standard Deviation should be equal or greater than 0');
                return
            end
            this.cam_sd = value;
            this.Cam.sd = value;
        end
        function set.cube_sd( this, value )
            if value < 0
                warning('cube_sd must be equal or greater than 0');
                return;
            end
            this.cube_sd = value;
            this.Cube = CHalfCube( 1, eye(3), zeros(3,1), this.cube_sd );
        end            
        
        function clearFigures( this )
            cla(this.hF_scene);
            cla(this.hF_image);
        end
        
        function updatePose( this )
            % Create new pose
            [R_w_c, t_w_c] = this.poseFactory.gen( 1 );
            this.R_gt = R_w_c';
            
            % Assign new camera pose
            this.Cam.R = R_w_c;
            this.Cam.t = t_w_c;
            
            if this.WITH_PLOT
                subplot(this.hF_scene);
                this.Cube.plotScene(this.Cam);
            end
        end
        
        function [err_elqursh, global_err_elqursh,...
                  err_OP3A, global_err_OP3A,...
                  global_err_OP3A_W] = ...
                  compareSingleCubeResults( this )
            [err_elqursh, global_err_elqursh] = this.syntheticElqursh;
            [err_OP3A, global_err_OP3A, global_err_OP3A_W] = this.syntheticOP3A;
            
            % Plot results
            h_ = figure;
            title('Error histogram'); hold on;
            
            subplot(1,2,1); hold on;
            plot(1, err_elqursh, '*b');
            plot(2, err_OP3A, '*g');
            ax = axis; axis([0 3 ax(3:4)]);
            
            subplot(1,2,2); hold on;
            hist( err_elqursh, 20, 0:0.1:10 );
            hist( err_OP3A, 0:0.1:10 );
            h = findobj(gca,'Type','patch');
            set(h(2),'FaceColor','b','EdgeColor','k');
            set(h(1),'FaceColor','g','EdgeColor','k');
            
            set(h_,'units','normalized','outerposition',[0 0 1 1]);
        end
        
        function [err_elqursh, global_err_elqursh,...
                  err_OP3A, global_err_OP3A,...
                  global_err_OP3A_W] = ...
                  compareSeveralCubeResults( this, Nsim )
            err_elqursh = cell(Nsim,1);
            err_OP3A = cell(Nsim,1);
            global_err_elqursh = zeros(1,Nsim);
            global_err_OP3A = zeros(1,Nsim);
            global_err_OP3A_W = zeros(1,Nsim);
            for ii = 1:Nsim
                this.updatePose;
                
                [err_elqursh{ii}, global_err_elqursh(ii)] = ...
                    this.syntheticElqursh;
                [err_OP3A{ii}, global_err_OP3A(ii), global_err_OP3A_W(ii)] = ...
                    this.syntheticOP3A;
            end
            % Linealize all values
            err_elqursh = [err_elqursh{:}];
            err_OP3A = [err_OP3A{:}];
            
            if nargout == 0
                % Store simulation data
                file = datestr(now);
                file = strrep( file,'-','_' );
                file = strrep( file,' ','_' );
                file = strrep( file,':','_' );
                file = ['Several_',file,'.mat'];
                file = fullfile(pwd,'Odometry','store',file);
                save(file, 'err_elqursh','global_err_elqursh',...
                    'err_OP3A','global_err_OP3A','global_err_OP3A_W' );
                
                this.plotSeveralCubeResults(...
                    err_elqursh, global_err_elqursh,...
                    err_OP3A, global_err_OP3A, global_err_OP3A_W );
            end
        end
        
        function plotSeveralCubeResults( ~,...
                  err_elqursh, global_err_elqursh,...
                  err_OP3A, global_err_OP3A, global_err_OP3A_W )
            
            % Plot results
            h_ = figure;
            fm = 1; fn = 3;
            title('Error histogram'); hold on;
            
            subplot(fm,fn,1); hold on;
            plot(1, err_elqursh, '*b');
            plot(2, err_OP3A, '*g');
            ax = axis; axis([0 3 ax(3:4)]);
            
            subplot(fm,fn,2); hold on;
            hist( err_elqursh, 20, 0:0.1:10 );
            hist( err_OP3A, 0:0.1:10 );
            h = findobj(gca,'Type','patch');
            set(h(2),'FaceColor','b','EdgeColor','k');
            set(h(1),'FaceColor','g','EdgeColor','k');
            
            subplot(fm,fn,3); hold on;
            hist( global_err_elqursh, 20, 0:0.1:10 );
            hist( global_err_OP3A, 0:0.1:10 );
            hist( global_err_OP3A_W, 0:0.1:10 );
            h = findobj(gca,'Type','patch');
            set(h(3),'FaceColor','b','EdgeColor','k');
            set(h(2),'FaceColor','g','EdgeColor','k');
            set(h(1),'FaceColor','r','EdgeColor','k');
            
            set(h_,'units','normalized','outerposition',[0 0 1 1]);
        end
        
        function testDistance( this, arr_d, Nsim )
            N_methods = 3;
            
            Nd = numel(arr_d);
            err1 = cell(1,Nd);
            err2 = cell(1,Nd);
            gerr1 = cell(1,Nd);
            gerr2 = cell(1,Nd);
            gerr3 = cell(1,Nd);
            
            for ii=1:Nd
                d = arr_d(ii);
                this.poseFactory.min_d = d;
                this.poseFactory.max_d = d;
                [err1{ii}, gerr1{ii}, err2{ii}, gerr2{ii}, gerr3{ii}] = ...
                    this.compareSeveralCubeResults( Nsim );
                
                Cval{ii} = num2str(d);
            end
            
            % Store simulation data
            file = datestr(now);
            file = strrep( file,'-','_' );
            file = strrep( file,' ','_' );
            file = strrep( file,':','_' );
            file = ['testDistance_',file,'.mat'];
            file = fullfile(pwd,'Odometry','store',file);
            save( file, 'err1','gerr1', 'err2','gerr2', 'gerr3',...
                'arr_d', 'Nsim' );
            aux_struct = struct( this ); %#ok<NASGU>
            save( file, '-struct', 'aux_struct',...
                    'cam_sd', 'cube_sd', 'Cam', 'Cube', '-append' );
            
            % Matrix of stacked data
%             M = [ cell2mat( gerr1' )', cell2mat( gerr2' )' ];
            M = [ cell2mat( gerr1' )', cell2mat( gerr2' )', ...
                  cell2mat( gerr3' )' ];
            
            Nx = Nd;
            % Cval contains the string with X-value
            % for each column of data
            Cval = repmat( Cval, 1, N_methods );
            % Ctag contains the string with method corresponding
            % to each column of data
%             Ctag = [ repmat({'Elqursh'},1,Nx),...
%                      repmat({'OP3A'},1,Nx) ];
            Ctag = [ repmat({'Elqursh'},1,Nx),...
                     repmat({'OP3A'},1,Nx),...
                     repmat({'OP3A-W'},1,Nx)];
            
            % Parameters to control the position in X label
            Npos    = 5;    % gap between samples in X label
            pos_ini = 1;    % initial value in X label
            Nsep    = 0.5;  % gap between methods in X label
            % Load the vector of positions
            pos_aux = pos_ini:Npos:Npos*Nx;
            pos_    = pos_aux;
            pos_1   = [];
            for i = 1:N_methods-1
                pos_ = [pos_ pos_aux+i*Nsep];
            end
            
            color_ = [0.2980392156862745 0.4470588235294118 0.6901960784313725;
                0.3333333333333333 0.6588235294117647 0.40784313725490196;
%                 0.7686274509803922 0.3058823529411765 0.3215686274509804;
                %0.5058823529411764 0.4470588235294118 0.6980392156862745;];
                0.8                0.7254901960784313 0.4549019607843137];
            color = repmat(color_,Nx,1);
            
            h = figure; hold on;
            boxplot(M,{Cval,Ctag},...
                'position',sort(pos_),'colors', color,...
                'factorgap',0,'whisker',0,'plotstyle','compact');
            
            % Remove the outliers
            bp_ = findobj(h, 'tag', 'Outliers'); % Find handler
            set(bp_,'Visible','Off'); % Remove object
            
            % Plot lines
            median_ = median(M);
            for i = 1:N_methods
                x_ = pos_(1, Nx*(i-1)+1:Nx*i);
                y_ = median_(1, Nx*(i-1)+1:Nx*i);
                plot(x_,y_,'Color',color(i,:),'LineWidth',1.5);
            end
            
            % Plot legend
            Cleg = {'Elqursh','O3PA','O3PA-W'};
            Clab = {Cval{1,1:Nx}};
%             set(gca,'YScale','log');
            set(gca,'XTickLabel',{' '});
            [legh,objh,outh,outm] = legend(Cleg);
            set(objh,'linewidth',3);
            set(gca,'XTick',pos_aux);
            set(gca,'XTickLabel',Clab);
        end
        
        function [err, global_err] = syntheticElqursh( this )
            
            % Assign variables from properties
            Cube = this.Cube;
            Cam  = this.Cam;
            R_gt = this.R_gt;
            
            if this.WITH_PLOT
                figure(this.hFigure);
            end
            
            % Solve with Elqursh
            primitives = Cube.elqursh_primitives;
            Nprim = length(primitives);
            K = Cam.K;
            err = zeros(1,length(primitives));
            I = eye(3);
            cell_d = cell(1,Nprim);
            cell_v = cell(1,Nprim);

            for k=1:Nprim
                prim = primitives{k};
                
                if this.WITH_PLOT
                    subplot(this.hF_scene);
                    h_prim = Cube.plot_prim( prim );
                end
                
                prim2D = prim.project(Cam);
                if this.WITH_PLOT
                    subplot(this.hF_image);
                    prim2D.plot;
                end
                
                line = cell(1,3);
                for ii=1:3
                    line{ii} = prim(ii).projectLine( Cam );
                end
                v2 = cross( line{2}, line{3} );
                v1 = null( [ v2'*(inv(K)'*inv(K)) ; line{1}' ] );
                
                % Approximate solution, only valid for almost canonical
                % directions
                [~,dir1] = max(abs(prim(1).v));
                [~,dir2] = max(abs(prim(2).v));
                dir3 = setdiff(1:3,[dir1 dir2]);
                
                R = eye(3);
                R(:,dir1) = snormalize(K \ v1);
                R(:,dir2) = snormalize(K \ v2);
                switch dir3
                    case 1
                        R(:,1) = cross(R(:,2),R(:,3));
                    case 2
                        R(:,2) = cross(R(:,3),R(:,1));
                    case 3
                        R(:,3) = cross(R(:,1),R(:,2));
                end
                
                % Compute 4 possible sign combinations
                CR = cell(1,4);
                CR{1} = R;
                CR{2} = R * diag([-1 1 -1]);
                CR{3} = R * diag([-1 -1 1]);
                CR{4} = R * diag([1 -1 -1]);
                % Find good one minimizing Frobenius norm
                d = zeros(1,4);
                for iR = 1:4
                    d(iR) = angularDistance(CR{iR},R_gt);
                end
                [~,iR] = min(d);
                R = CR{iR};
                % Error
                err(k) = angularDistance(R,R_gt);
                
                % Store vanishing points and directions for global Procrustes
                cell_d{k} = I(:,[dir1 dir2]);
                %     cell_d{k} = [prim(1).v, prim(2).v];
                cell_v{k} = R(:,[dir1 dir2]); % Equivalent to inv(K)*K*rk
            end
            % Set camera image borders
            if this.WITH_PLOT
                subplot(this.hF_image);
                Cam.setImageBorder;
            end
            
            D1 = [cell_d{:}];
            D2 = [cell_v{:}];
            
            % Compute through Procrustes problem:
            [U,~,V] = svd( D2*D1' );
            R_pro = U*V';
            global_err = angularDistance(R_pro,R_gt);
        end
        
        function [err, global_err, global_err_W] = syntheticOP3A( this )
            
            Cube = this.Cube;
            Cam  = this.Cam;
            R_gt = this.R_gt;
            
            if this.WITH_PLOT
                figure(this.hFigure);
            end
            
            % Compute with Trihedron method
            primitives = Cube.trihedron_primitives;
            Nprim = length(primitives);
            
            K = Cam.K;
            err = zeros(1,Nprim);
            cell_d = cell(1,Nprim);
            cell_v = cell(1,Nprim);

            for k=1:Nprim
                prim = primitives{k};
                % Take care of segments orientation
                % (building a dextrorotatory system)
                for kk=1:3
                    % Valid only for canonical directions (simulation)
                    if any(prim(kk).v < 0)
                        prim(kk) = prim(kk).inverse;
                    end
                end
                
                if this.WITH_PLOT
                    subplot(this.hF_scene);
                    h_prim = Cube.plot_prim( prim );
                end
                
                prim2D = prim.project(Cam);
                if this.WITH_PLOT
                    subplot(this.hF_image);
                    prim2D.plot('-k',{'x','y','z'});
                end
                
                for ii=1:3
                    l = snormalize( prim(ii).projectLine( Cam ) ); % Image coordinates of line
                    
                    cell_n{ii} = snormalize( K' * l );
                    if this.WITH_COVARIANCE
                        objs_l(ii) = Manifold.S2( snormalize(l) );
                        
                        Am_l = this.lin_sd^2 * eye(2);
                        objs_l(ii).setMinimalCov( Am_l );
                        A_l = null(l') * Am_l * null(l')';
                        J_n_l = Dsnormalize(K'*l) * K';
                        objs_n(ii) = Manifold.S2( cell_n{ii} );
%                         objs_n(ii).setRepresentationCov( J_n_l * A_l * J_n_l' );
                        % Temporary
                        sd_n = 1e-6;
                        objs_n(ii).setMinimalCov( sd_n * eye(2) );
                    end
                end
                
                Nbp = [cell_n{:}];
                trihedronSolver = CTrihedronSolver( Nbp, K );
                trihedronSolver.loadSegments( prim2D );
                V_tri = trihedronSolver.solve;
                if this.WITH_COVARIANCE
                    trihedronSolver.loadCovariance( ...
                        {objs_n.A_X} );
                    trihedronSolver.computeCovariance;
                    obj_V = trihedronSolver.obj_V;
                    
                    WITH_MONTECARLO = false;
                    if WITH_MONTECARLO
                        keyboard
                        % Set manifold framework inputs
                        obj_L = Manifold.Dyn( objs_l(1), objs_l(2), objs_l(3) );
                        A_L = blkdiag( objs_l.A_x );
                        obj_L.setMinimalCov( A_L );
                        
                        obj_N = Manifold.Dyn( objs_n(1), objs_n(2), objs_n(3) );
                        obj_N.setMinimalCov( blkdiag( objs_n.A_x ) );
                        
                        % Set temporary necessary variables
                        this.prim2D = prim2D;

%                         out = Manifold.MonteCarloSim( ...
%                             @(X)this.Fun_OP3A_Fast(X),...
%                             obj_L, 'Ref', obj_V, 'N', 1e3 );
                        out = Manifold.MonteCarloSim( ...
                            @(X)this.Fun_OP3A_Fast_N(X),...
                            obj_N, 'Ref', obj_V, 'N', 1e4 );
                    end
                end
                
%                 if this.WITH_PLOT
%                     plotHomLineWin(K'\trihedronSolver.Nbp(:,3), 'r');
%                 end
                
                % Error
                err(k) = angularDistance(V_tri,R_gt);
                
                % Store vanishing points and directions for global Procrustes
                cell_d{k} = eye(3);
                cell_v{k} = V_tri; % Equivalent to inv(K)*K*rk
                
                % Store weight according to covariance
                cell_w{k} = trace( obj_V.A_x ) * eye(3);
            end
            % Set camera image borders
            if this.WITH_PLOT
                subplot(this.hF_image);
                Cam.setImageBorder;
            end
            
            D1 = [cell_d{:}];
            D2 = [cell_v{:}];
                       
            % norm(D1-R'*D2,'fro')
            % Compute through Procrustes problem:
            [U,~,V] = svd( D2*D1' );
            R_pro = U*V';
            global_err = angularDistance(R_pro,R_gt);
            
            % Weight matrix
            W  = blkdiag( cell_w{:} );
%             W  = W / trace(W);
            D1W = D1 / W;
            D2W = D2 / W;
            
            % Compute through weighted Procrustes problem:
            [U,~,V] = svd( D2W*D1W' );
            R_pro_W = U*V';
            global_err_W = angularDistance(R_pro_W,R_gt);
        end
        
        function V = Fun_OP3A_Fast( this, L )
            L = reshape( L.X, 3,3 );
            N = snormalize( this.Cam.K' * L );
            
            trihedronSolver = CTrihedronSolver( N, this.Cam.K );
            trihedronSolver.loadSegments( this.prim2D );
            V = trihedronSolver.solve;
            V = Manifold.SO3( V );
        end
        function V = Fun_OP3A_Fast_N( this, N )
            N = reshape( N.X, 3,3 );
            
            trihedronSolver = CTrihedronSolver( N, this.Cam.K );
            trihedronSolver.loadSegments( this.prim2D );
            V = trihedronSolver.solve;
            V = Manifold.SO3( V );
        end
    end
    
end

