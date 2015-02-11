% Script to generate experimentation
experiments = {'cam_sd', 'cube_sd', 'd', 'FOV', 'Rrel_sd', 'trel_sd'};
xlabels = {'\sigma_c (pixels)', '\sigma_L (m)', 'd (m)',...
           'FOV (deg)', '|R_{rel}| (deg)', '|t_{rel}| (m)'};
xvec = {logspace(-3,0,10),... %[0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1],...
        logspace(-3,-1,10),... %[0.001, 0.005, 0.01, 0.05],...
        3:2:17,...
        700:-100:200,...
        [0.1, 1, 2, 5, 10, 20],...
        [0.001, 0.01, 0.1, 1]};
Nsim = 1000;

figures_path = '/home/jesus/P3oA_2014/figures/experiments';
for i=1:numel(experiments)
    fprintf('Launching %s experiment\n',experiments{i});
    test = COdometryTest(false); test.WITH_DEBUG = false;
    if strcmp(experiments{i},'cube_sd')
        test.Cam.sd = 0; % Standard noise to apply in cube case
    else
        test.Cam.sd = 0.1; % Standard noise to apply in all cases
    end
    
    test.test( @(x)test.update(experiments{i},x), xvec{i}, Nsim );
    xlabel(xlabels{i});
    ylabel('Rotation error (deg)');
%     title(sprintf('Experiment %s',experiments{i}));
    figname = ['Synth_Experiment_',experiments{i}];
    set(gca,'YScale','log');
    print(gcf, '-depsc', fullfile(figures_path,figname));
    clear test;
end