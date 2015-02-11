function plotComparisonChart( comp_obj, cam_sd_vals, scan_sd_vals, Nobs_vals, with_outliers )
% Plot options
plot_sim_file = fullfile( pwd, 'plotBoxplot.ini' );
plotOpts = readConfigFile( plot_sim_file );
extractStructFields( plotOpts );
clear plotOpts;

% Check size of inputs (only one should be > 1)
if ~exist( 'cam_sd_vals', 'var' )
    cam_sd_vals = comp_obj.cam_sd_vals;
end
if ~exist( 'scan_sd_vals', 'var' )
    scan_sd_vals = comp_obj.scan_sd_vals;
end
if ~exist( 'Nobs_vals', 'var' )
    Nobs_vals = comp_obj.Nobs_vals;
end
if ~exist( 'with_outliers', 'var' )
    with_outliers = true;
end

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
cam_idxs  = comp_obj.getIndexes('cam_sd_vals',cam_sd_vals);
scan_idxs = comp_obj.getIndexes('scan_sd_vals',scan_sd_vals);
Nobs_idxs = comp_obj.getIndexes('Nobs_vals',Nobs_vals);

% Extract dim x Nsim matrices for representation
all_R_err = [];
all_t_err = [];
Cval = cell(1,0);
Ctag = cell(1,0);
Cleg = cell(1,0);
Clab = cell(1,0);
N_met = 0;
for idx_pattern = 1:size(patterns,1)
    pat = patterns{idx_pattern,1};
    methods = patterns{idx_pattern,2};
    for idx_method = 1:length(methods)
        met = methods{idx_method};
        %                     Rt = squeeze(obj.(pat).(met).mem(:,scan_sd_it,N_co_it,:));
        Rt = squeeze(comp_obj.(pat).(met).mem(cam_idxs,scan_idxs,Nobs_idxs,:));
        if size(Rt,2)==1
%             warning('Check squeeze changed dimension order if too many singletons');
            Rt = Rt';
        end
        R_err = cell2mat( cellfun(@(x)comp_obj.R_err(x), Rt, 'UniformOutput',false) );
        t_err = cell2mat( cellfun(@(x)comp_obj.t_err(x), Rt, 'UniformOutput',false) );
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
err.xtick = num2str( comp_obj.cam_sd_vals );

% Parameters to control the position in X label
Npos    = 5;    % gap between samples in X label
pos_ini = 1;    % initial value in X label
Nsep    = 0.5;  % gap between methods in X label
% Load the vector of positions
pos_aux = pos_ini:Npos:Npos*Nx;
pos_    = pos_aux
pos_1   = [];
for i = 1:N_met-1
    pos_ = [pos_ pos_aux+i*Nsep];
end

color_ = [0.2980392156862745 0.4470588235294118 0.6901960784313725;
    0.3333333333333333 0.6588235294117647 0.40784313725490196;
    0.7686274509803922 0.3058823529411765 0.3215686274509804;
    %0.5058823529411764 0.4470588235294118 0.6980392156862745;];
    0.8                0.7254901960784313 0.4549019607843137];
color = repmat(color_,Nx,1);
%             color = repmat(rand(N_met,3),Nx,1);
Rlab = 'Rotation error (deg)';
tlab = 'Translation error (m)';
Cleg = {'Trihedron','Kwak et al.','Wasielewski et al.','Vasconcelos et al.'};

% Boxplot for R
h = figure; hold on;
boxplot(all_R_err,{Cval,Ctag},'position',sort(pos_),'colors', color, 'factorgap',0,'whisker',0,'plotstyle','compact');
if ~with_outliers
    bp_ = findobj(h, 'tag', 'Outliers');
    set(bp_,'Visible','Off');   % Remove the outliers
end
xlabel(xlab); ylabel(Rlab);
% Plot the lines
median_ = median(all_R_err,1);
for i = 1:N_met
    x_ = pos_(1, Nx*(i-1)+1:Nx*i);
    y_ = median_(1, Nx*(i-1)+1:Nx*i);
    plot(x_,y_,'Color',color(i,:),'LineWidth',1.5);
    %Cleg = {Cleg{:}, Ctag{1,Nx*(i-1)+1} };
end
Clab = {Cval{1,1:Nx}};
set(gca,'YScale','log');
set(gca,'XTickLabel',{' '});
[legh,objh,outh,outm] = legend(Cleg);
set(objh,'linewidth',3);
set(gca,'XTick',pos_aux);
set(gca,'XTickLabel',Clab);

% Boxplot for t
h = figure; hold on;
boxplot(all_t_err,{Cval,Ctag},'position',sort(pos_),'colors', color, 'factorgap',0,'whisker',0,'plotstyle','compact');
if ~with_outliers
    bp_ = findobj(h, 'tag', 'Outliers');
    set(bp_,'Visible','Off');   % Remove the outliers
end
xlabel(xlab); ylabel(tlab);
% Plot the lines
median_ = median(all_t_err, 1);
for i = 1:N_met
    x_ = pos_(1, Nx*(i-1)+1:Nx*i);
    y_ = median_(1, Nx*(i-1)+1:Nx*i);
    plot(x_,y_,'Color',color(i,:),'LineWidth',1.5);
end
set(gca,'YScale','log');
set(gca,'XTickLabel',{' '});
[legh,objh,outh,outm] = legend(Cleg);
set(objh,'linewidth',3);
set(gca,'XTick',pos_aux);
set(gca,'XTickLabel',Clab);
end