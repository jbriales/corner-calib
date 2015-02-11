function [err_R] =...
    testOrientation( this, dataset_folder, dataset_type,...
    i0, step, debug )
if ~exist('i0','var')
    i0 = 1;
end
if numel(i0)==1
    i0 = [i0 +inf];
end
if ~exist('step','var')
    step = 1;
end

if ~exist('debug','var')
    debug = false;
end

raw = CRawlogCam( dataset_folder, [], dataset_type );

% Get camera configuration
config_file = fullfile( pwd, 'configs', [raw.typeOfSource,'.ini'] );
SConfigCam = readConfigFile( config_file );
this.Cam = CRealCamera( SConfigCam );

tracker{1} = CGtracker([],false);
tracker{2} = CGtracker([],false);

% Set GT triplets (for P3oA method) (by visual inspection)
switch raw.typeOfSource
    case 'Freiburg'
        % Freiburg
        idxs_GT = {[1 4 8 10]
            [2 5 7 12]
            [3 6 9 11]};
    case 'Blender'
        % Blender
        idxs_GT = {[1 4 7]
            [2 5 8]
            [3 6 9]};
        % Get 3D GT segments
        cube_GT = CHalfCube( 2, eye(3), ones(3,1) );
        segs3D_GT = [cube_GT.segments{:}];
        CGT.segs3D( segs3D_GT );
        %                     figure; hold on;
        %                     segs3D_GT.plot3('-k',{'1','2','3 ','4','5','6','7','8','9'});
        %                     axis equal
        % Conversion map
        
end
triplets_GT = allcomb( idxs_GT{:} )';
% Use fixed ordering: ascend order (since x-y-z id does not matter)
triplets_GT = sort(triplets_GT,1,'ascend'); %#ok<UDIM>
% Order columns
triplets_GT = sortrows(triplets_GT')';
% DEBUG: Store persistent value
CGT.triplets(triplets_GT);

% Create combinations of all frames in chosen interval
if 0
frame_pairs = nchoosek(i0(1):step:i0(2),2)';
else
frame_pairs = nchoosek(i0(1):i0(2),2)'; % Take all combinations
to_remove = abs(frame_pairs(1,:)-frame_pairs(2,:))<step;
frame_pairs(:, to_remove) = [];
end
N_frame_pairs = size(frame_pairs,2);

% Treat err as a struct
err_R = repmat( struct('P3oA',[],'Elqu',[],'Fuse',[]), 1, N_frame_pairs );
vec_t_norm = NaN( 1, N_frame_pairs );
vec_R_norm = NaN( 1, N_frame_pairs );

% Preallocate image handles
hSegs = cell(1,2);
gSegs = cell(1,2);

iptsetpref('ImshowBorder', 'tight');

for i=1:N_frame_pairs
    fr = zeros(1,2);
    fr(1) = frame_pairs(1,i);
    fr(2) = frame_pairs(2,i);
    
    for k=1:2
        % Load segments for current step
        tracker{k}.loadSegs( raw.frames(fr(k)).path_metafile );
        
        if debug
            % Load image
            tracker{k}.loadImage( raw.frames(fr(k)).loadImg );
        end
    end
    
    if debug
    % Show both frames 
    if ~exist('hIm','var')
        % First time
        figure('Name','Track figure');
%         hIm = imshowpair( tracker{1}.img, tracker{2}.img, 'montage' ); hold on;
        hIm = imshowpair( tracker{1}.img, tracker{2}.img, 'falsecolor' ); hold on;
        freezeColors;
        %                     set(gcf,'Position',...)
    else
        % Update color data (use imfuse for getting composite image)
        set(hIm,'CData',...
            imfuse( tracker{1}.img, tracker{2}.img, 'falsecolor' ) );
    end
    
    % Show segments and tags
    for k=1:2
        tags{k} = tracker{k}.segs.tags( tracker{k}.maskSegs );
%         if exist('hSegs','var')
        if ishandle(hSegs{k})
            delete(hSegs{k});
            delete(gSegs{k});
        end
        if 0
            hSegs{k} = tracker{k}.segs.segs(tracker{k}.maskSegs).plot('r');
            % Meeting case: 7-8-9
            set(hSegs{k}([7 8 9]),'Color',[255 128 0]/255,'LineWidth',2);
            % Non-meeting case: 1-3-5
            set(hSegs{k}([1 3 5]),'Color',[0 0 204]/255,'LineWidth',2);
        else
        [hSegs{k}, gSegs{k}] = tracker{k}.segs.segs(tracker{k}.maskSegs).plot('r', tags{k});
        % Change tags color for better visualization
        set(gSegs{k},'Color','r');
        end
    end
    end
    
    % Extract and process data
    for k=1:2
        % Extract segments arrays (to avoid accessing problems)
        segs{k} = tracker{k}.segs.segs;
    end
    [common_tags, I1, I2] = intersect( ...
        [segs{1}(tracker{1}.maskSegs).tag],...
        [segs{2}(tracker{2}.maskSegs).tag] ); %#ok<ASGLU>
    % Create matches array: 2xM matrix where each column
    % [idx1;idx2] gives the indexes of segments in arrays
    % segs{1} and segs{2} that form a potential match pair
    %                 matches = [I1, I2]'; % For general cases
    matches = repmat( common_tags, 2,1 ); % For already matched case
    
    % Get GT relative pose from Freiburg
    R_gt = raw.frames(fr(1)).pose.R' * raw.frames(fr(2)).pose.R;
    CGT.R(R_gt); % Store persistent GT value
    t_gt_norm = norm(raw.frames(fr(2)).pose.t - raw.frames(fr(1)).pose.t );
                
    %% Compute relative values
    % Testing with Elqursh code
    [R_Elqu, V_Elqu] = this.codeElqursh( segs{1}, segs{2}, matches',...
        this.Cam.size, tracker{1}, tracker{2} );
    R_Elqu = R_Elqu'; % Transpose because according to Elqursh notation R = (^c2)R(_c1)
    
    % Testing with P3oA code
    % Rotation threshold depends a lot on dataset error (higher
    % in reality than Blender)
    switch raw.typeOfSource
        case 'Freiburg'
%             rotThres = 3; % For old RANSAC version
            rotThres = 1.5; % For new RANSAC version
        case 'Blender'
            rotThres = 0.5;
    end
    % DEBUG:
    if exist('cube_GT','var')
        CGT.any('Cam1',Cam1);
        CGT.any('Cam2',Cam2);
        CGT.any('Cube',cube_GT);
    end
    [R_P3oA, V_P3oA, triplets] = this.codeP3oA( segs{1}, segs{2}, matches, this.Cam.K, rotThres );
%     fprintf('Compare GT and estimate rotation [deg]: %f\n',angularDistance(CGT.R,R_P3oA) );
    
    % Fuse data from both methods
    V1 = [V_P3oA{1}, V_Elqu{1}];
    V2 = [V_P3oA{2}, V_Elqu{2}];
    % Test for outliers in fused set
    d = acosd( dot( V1, R_P3oA*V2, 1) );
    %                 d = asind( sqrt( sum( (V1 - R_P3oA*V2).^2, 1 ) ) );
    %                 inliers = find( d < rotThres );
    inliers = d < rotThres;
    if this.WITH_DEBUG
        fprintf('Inliers for Fused: %d/%d\n',sum(inliers),size(V1,2));
    end
    % Keep only inliers
    V1 = V1(:,inliers);
    V2 = V2(:,inliers);
    
    [Ur,~,Vr] = svd(V1*V2');
    if (det(Ur)*det(Vr)>=0), Sr = eye(3);
    else Sr = diag([1 1 -1]);
    end
    R_Fuse = Ur*Sr*Vr';
    
    %% Compute errors
    % Store translation magnitude
    vec_t_norm(i) = t_gt_norm;
    % Store rotation magnitude
    vec_R_norm(i) = angularDistance(R_gt,eye(3));
    % Relative errors in rotation
    err_R(i).P3oA = angularDistance(R_P3oA,R_gt);
    err_R(i).Elqu = angularDistance(R_Elqu,R_gt);
    err_R(i).Fuse = angularDistance(R_Fuse,R_gt);
    
    % Print (display) values
    fprintf('#It\tFrames\t|R_GT|[deg]\t|t|[cm]\tP3oA\tElqursh\tFuse\n');
    fprintf('%d/%d\t%d-%d\t', i, N_frame_pairs, fr(1), fr(2));
    fprintf('%f\t%.2f\t%.4f\t%.4f\t%.4f\n',...
            angularDistance(R_gt,eye(3)),t_gt_norm*100,...
            err_R(i).P3oA, err_R(i).Elqu, err_R(i).Fuse );
    if ~isempty( CGT.missedInliers( triplets ) )
    fprintf('P3oA rotation missed inliers:\n');
    disp( CGT.missedInliers( triplets ) );
    end
    if ~isempty( CGT.slippedOutliers( triplets ) )
    fprintf('P3oA rotation slipped outliers:\n');
    disp( CGT.slippedOutliers( triplets ) );
    end
end
fprintf('Finished test for frames #%d to #%d with step %d\n',i0(1),i0(2),step);

method = 'FreiRot';

% Store experiment data
file = datestr(now);
file = strrep( file,'-','_' );
file = strrep( file,' ','_' );
file = strrep( file,':','_' );
file = ['test_',method,'_',file,'.mat'];
file = fullfile(pwd,'Odometry','store',file);
save( file, 'err_R', 'method' );

% keyboard
%% Plot statistics of errors
% Colors configuration
c_blue = [51 51 255]/255;
c_oran = [255 128 0]/255;
c_gree = [102 255 102]/255;
c_face_blue = [153 153 255]/255;
c_face_oran = [255 204 153]/255;
c_face_gree = [102 255 102]/255;
color = {c_blue, c_oran};
color_face = {c_face_blue,c_face_oran};

figures_path = '/home/jesus/P3oA_2014/figures/experiments';

m = 1; n = 2;
figure('Name','Error of results'); hold on;
% set(gcf, 'renderer', 'opengl')
% Boxplots
% subplot(m,n,1); hold on;
% title('Boxplot of rotation angle errors');
M_R = [ [err_R.P3oA]
        [err_R.Elqu] ]';
boxplot(M_R, {'P3oA','Elqursh'},...
        'boxstyle','outline',...
        'whisker',1.5,...
        'factorgap',0,...
        'positions',[0.5 1]);

% Style parameters
LW = 2; % Line Width
hBox = sort( findobj(gca,'Tag','Box') );
hOut = sort( findobj(gca,'Tag','Outliers') );
for k=1:2
    hFill(k) = patch(get(hBox(k),'XData'),get(hBox(k),'YData'),'k');
    set(hFill(k),'FaceColor',color_face{k});
%     set(hFill(k),'FaceAlpha',0.5); % Only with OpenGL renderer
    set(hFill(k),'EdgeColor',color{k});
    set(hFill(k),'LineWidth',LW);
    
    set(hOut(k),'Marker','.');
    set(hOut(k),'Color',color{k});
end
Tags = {'Median','Lower Adjacent Value','Upper Adjacent Value',...
        'Lower Whisker','Upper Whisker'};
for i=1:numel(Tags)
    h_ = sort( findobj(gca,'Tag',Tags{i}) );
    for k=1:2
        set(h_(k),'Color',color{k},...
                  'LineWidth',LW,...
                  'LineStyle','-');
    end
end
hMed = sort( findobj(gca,'Tag','Median') );
for k=1:2
%     uistack(hMed(k),'bottom');
    h__ = plot( get(hMed(k),'XData'), get(hMed(k),'YData') );
    set(h__,'Color',color{k},...
            'LineWidth',LW,...
            'LineStyle','-');
end
ylabel('Rotation error (deg)');

figname = ['Real_Experiment_',num2str(step),'_box'];
print(gcf, '-depsc', fullfile(figures_path,figname));

%     'Outliers'
%     'Median'
%     'Box'
%     'Lower Adjacent Value'
%     'Upper Adjacent Value'
%     'Lower Whisker'
%     'Upper Whisker'

% Histograms
figure
% subplot(m,n,2); hold on;
% title('Histogram of rotation angle errors');
% Continuous approximation
% % ksdensity( M_R(:,1) )
% % ksdensity( M_R(:,2) )
% % h = findobj(gca,'Type','line');
% % set(h(2),'Color',c_blue); % Elqursh
% % set(h(1),'Color',c_oran); % P3oA
% Discrete histogram
xbinscenters = 0:0.1:2.5;
[N] = hist( [[err_R.P3oA]
             [err_R.Elqu] ]', xbinscenters );
% bar(xbinscenters,N./repmat(trapz(xbinscenters,N),size(N,1),1),'histc');
bar(xbinscenters,N,'histc');
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
axis_ = axis;
axis_(2) = 2.5;
axis(axis_);

if 0
figname = ['Real_Experiment_',num2str(step),'_box'];
print(gcf, '-depsc', fullfile(figures_path,figname));
end

if 0
% Ordered representation wrt |trel|
subplot(m,n,3); hold on;
title('Representation of errors wrt |t_{rel}|');
[dt,Idt] = sort(vec_t_norm,2,'ascend');
% [N,bin] = histc(dt, [0.01 0.1 0.2 0.5]);
plot( dt,[err_R(Idt).P3oA], 'r.' );
plot( dt,[err_R(Idt).Elqu], 'g.' );
% plot( dt,[err_R(Idt).Fuse], 'b.' );

% Ordered representation wrt |Rrel|
subplot(m,n,4); hold on;
title('Representation of errors wrt |R_{rel}|');
[dR,IdR] = sort(vec_R_norm,2,'ascend');
plot( dR,[err_R(IdR).P3oA], 'r.' );
plot( dR,[err_R(IdR).Elqu], 'g.' );
% plot( dR,[err_R(IdR).Fuse], 'b.' );
end

% keyboard
end