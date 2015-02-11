% Script to test numerical stability

% Add necessary toolboxes
addpath('/home/jesus/Dropbox/Computer_Vision/Matlab_tests/VisionBase');
addpath(genpath('/home/jesus/Dropbox/Computer_Vision/Matlab_tests/PeterKoversi'));
cd /home/jesus/Code/Corner_Calibration
addpath(genpath(pwd))

if 1
    test = COdometryTest(false); test.WITH_DEBUG = false;
    test.Cam.sd = 0;

    arr_x = 1;
    Nsim = 1e5;

    err_R = test.test( 'num', arr_x, Nsim );
else
test = COdometryTest(false); test.WITH_DEBUG = false;
test.Cam.sd = 0;

Nsim = 100;

err_R = repmat( struct('P3oA',[],'Elqu',[]), Nsim, 1 );
status = 0;
counter= 0;
total  = Nsim;
tic
for ii = 1:Nsim
    test.updatePose; % Update in each iteration
    
    try
        err_R(ii) = test.syntheticOrientation;
    catch exception
        disp(exception.message);
        for kstack=1:numel(exception.stack)
            disp(exception.stack(kstack));
        end
        fields = fieldnames(err_R);
        for k=1:numel(fields)
            f = fields{k};
            err_R(ii).(f) = NaN;
        end
    end
    % Print status
    counter = counter + 1;
    if floor(counter/total*100) > status
        status = floor(counter/total*100);
        fprintf('%d/100\t%d of %d\tEstimated remaining %d\n',...
            status,counter,total,floor((total-counter)*toc/counter));
    end
end

method = 'num';
% Store simulation data
file = datestr(now);
file = strrep( file,'-','_' );
file = strrep( file,' ','_' );
file = strrep( file,':','_' );
file = ['test_',method,'_',file,'.mat'];
file = fullfile(pwd,'Odometry','store',file);
save( file, 'err_R', 'method' );

keyboard
end

% Histograms
c_blue = [51 51 255]/255;
c_oran = [255 128 0]/255;
c_gree = [102 255 102]/255;
c_face_blue = [153 153 255]/255;
c_face_oran = [255 204 153]/255;
c_face_gree = [102 255 102]/255;
color = {c_blue, c_oran};
color_face = {c_face_blue,c_face_oran};
figure

M = [[err_R.P3oA]
     [err_R.Elqu] ]'; % Data as columns

edges = logspace(-15,-11,5);
% edges = logspace(-11,0,12);
N = histc(M,edges);
bar(N);
%# fix the x-labels, x-axis extents
centers = sqrt(edges(1:end-1).*edges(2:end));
xlim([0.5,length(centers)+0.5])
set(gca,'xticklabel',num2str(centers(:),'%5.2E'))
set(gca,'XTickLabel',...
        cellfun(@(x)['10^{',num2str(x),'}'],...
                num2cell(log10(edges(1:end-1))),'UniformOutput',false));
            
            
xbinscenters = logspace(-6,0);
[N,X] = hist( [[err_R.P3oA]
               [err_R.Elqu] ]', xbinscenters );
% bar(xbinscenters,N./repmat(trapz(xbinscenters,N),size(N,1),1),'histc');
bar(X,N,'histc');
h = sort(findobj(gca,'Type','patch'));
for k=1:2
    set(h(k),'FaceColor',color{k},'EdgeColor','w');
end
% set(h(2),'FaceColor',c_oran,'EdgeColor','w'); % P3oA
% set(h(1),'FaceColor',c_gree,'EdgeColor','w'); % Elqu
% set(h(1),'FaceColor',c_gree,'EdgeColor','w'); % W-P3oA
legend('P3oA','Elqursh','Location','NorthEast')
xlabel('Rotation error (deg)');
ylabel('Occurrences');
% axis_ = axis;
% axis_(2) = 2.5;
% axis(axis_);