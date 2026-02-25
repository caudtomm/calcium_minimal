
dbstop if error

s = true; % save figures to files?
savepath = 'bin1';
saveType = 'vector'; % 'vector' or 'raster'

% for sliding windows
window_duration = .3; % [seconds]
t_lim_sec = [-5 35]; % from 5 sec before to 35 seconds after stimulus onset
overlap = .1; % [seconds]

%% Load dataset
filename = 'odorexp004_IC1_130625.mat';

%
experiment = load(filename).a; % Experiment object
% or
% experiment = a; clear a

% avoid any spelling mismatches
for i = 1:numel(experiment.traces)
experiment.traces{i}.subject_group = experiment.subjectTab.group{i};
end


%% initialize output figure saving
cfg = PlotConfig('colormapName','lapaz','favouriteColors',[84,85,73]); % (test1, test2, ctrl)
cfg.savePath = savepath;

v = ExperimentViewer(experiment);
v.plotConfig = cfg;
cfg = v.plotConfig;

v.dataFilter.traceType = 'pSpike';
v.dataFilter.interval = [1,20];
v.dataFilter.trial_sorting = 'stim_id';
dft = v.dataFilter;
 

%% FIGURE 1

%% topography of the activity

% different single-unit metrics vs anatomy
out = v.plotTopographicHead('plotType','maps','method','avg intensity')
out = v.plotTopographicHead('plotType','maps','method','variance')
out = v.plotTopographicHead('plotType','maps','method','selectivity of tuning')
out = v.plotTopographicHead('plotType','maps','method','stability of tuning')
out = v.plotTopographicHead('plotType','maps','method','general suppression score')


% response cosine sim vs anatomy
v.dataFilter.interval = [.5 20]; % stimulus
[~,events] = ModeSelector(v).extract;
res = {};
for n = 1:42
    ta = TopographicAnalysis(v.filtered_traces{n});
    data = ActivityTraces.format(events{n});
    res{n} = ta.quantify(data,'cosine',100,10);
    title(['fish #',num2str(n)])
end
mantelRs = cellfun(@(x) x.mantelR,res);
figure; histogram(mantelRs,10)

% noise cosine sim vs anatomy
v.dataFilter.interval = [-22 -2]; % pre-stimulus
[~,events] = ModeSelector(v).extract;
res = {};
for n = 1:42
    ta = TopographicAnalysis(v.filtered_traces{n});
    data = ActivityTraces.format(events{n});
    res{n} = ta.quantify(data,'cosine',100,10);
    title(['fish #',num2str(n)])
end
mantelRs = cellfun(@(x) x.mantelR,res);
figure; histogram(mantelRs,10)

%% naive intertrial baseline correlation
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'naïve';
v.dataFilter.stims_allowed = 'all stimuli';
v.dataFilter.interval = [-22 -2];
v.dataFilter.trial_sorting = 'chronological';
hf = figure;
subplot(121)
C = v.plotDistancesHead;
C= 1-C.distMat3d;
[~,ntrials,nsubjects] = size(C);
y = [];
for i = 1:ntrials-1
    idx = triu(true(ntrials),i)-triu(true(ntrials),i+1);
    idx = logical(repmat(idx,1,1,nsubjects));
    y_avg = mean(C(idx),'omitmissing');
    y_std = std(C(idx),[],'omitmissing');
    y = [y; y_avg y_std];
end
subplot(122)
plotLineNShade(1:ntrials-1,y(:,1),y(:,2),'k',cfg)
ylim([0 1])
xlabel('Distance (trials)'); ylabel('r')
cfg.figSize = 'large';
cfg.aspRatioType = 'wide';
cfg.setFigure;
cfg.saveFigure(gcf,'naive base intertrial corr', saveType)
cfg = v.plotConfig;
v.dataFilter = dft;

%% corr between stimulus activity and post - prestimulus activity
post_interval = [115.5 135];
pre_interval = [-22.5 -2];
stim_interval = [.5 20];
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'trained';
v.dataFilter.stims_allowed = 'all stimuli';
v.dataFilter.trial_sorting = 'chronological';

% post vs pre
v.dataFilter.interval = post_interval;
[~,post_activity,all_labs] = ModeSelector(v).extract;
labs = all_labs{1};
post_activity = cellfun(@(x) squeeze(mean(x,1,'omitmissing')),post_activity,'UniformOutput',false);
v.dataFilter.interval = pre_interval;
[~,pre_activity] = ModeSelector(v).extract;
pre_activity = cellfun(@(x) squeeze(mean(x,1,'omitmissing')),pre_activity,'UniformOutput',false);
nsubjects = numel(pre_activity);
ntrials = width(pre_activity{1});
y = zeros(ntrials,nsubjects);
for i = 1:nsubjects
    for j = 1:ntrials
        prepattern = pre_activity{i}(:,j);
        postpattern = post_activity{i}(:,j);
        
        c = corrcoef(prepattern,postpattern);
        y(j,i) = c(2);
    end
end
t = 1:ntrials;
mu = mean(y,2,'omitmissing');
err = std(y,[],2,'omitmissing');
figure;
subplot(221)
plotLineNShade(t,mu,err,'k',cfg);
xticks(1:ntrials); xticklabels(labs)
axis tight; ylim([-1 1])
ylabel('post vs pre')
subplot(222)
histogram(y(:),50,'FaceColor','k','EdgeAlpha',0,'Orientation','horizontal')
hold on
line([0 max(xlim)],mean(y(:),'omitmissing')*[1 1],'Color','r','LineWidth',2)
ylabel('r')
ylim([-1 1])

% (post-pre) vs stim
basediff = cellfun(@(a,b) b-a, pre_activity, post_activity, 'UniformOutput',false);
v.dataFilter.interval = stim_interval;
[~,stim_activity] = ModeSelector(v).extract;
stim_activity = cellfun(@(x) squeeze(mean(x,1,'omitmissing')),stim_activity,'UniformOutput',false);
y = zeros(ntrials,nsubjects);
for i = 1:nsubjects
    for j = 1:ntrials
        diffpattern = basediff{i}(:,j);
        stimpattern = stim_activity{i}(:,j);
        
        c = 1-pdist([diffpattern,stimpattern]','correlation');
        y(j,i) = c;
    end
end
t = 1:ntrials;
mu = mean(y,2,'omitmissing');
err = std(y,[],2,'omitmissing');
subplot(223)
plotLineNShade(t,mu,err,'k',cfg);
xticks(1:ntrials); xticklabels(labs)
axis tight; ylim([-1 1])
ylabel('(post-pre) vs stim')
subplot(224)
histogram(y(:),50,'FaceColor','k','EdgeAlpha',0,'Orientation','horizontal')
hold on
line([0 max(xlim)],mean(y(:),'omitmissing')*[1 1],'Color','r','LineWidth',2)
ylabel('r')
ylim([-1 1])

set(gcf,'color','w')


%% naive corr with prestimulus in time
windows = defineTimeWindows(window_duration,t_lim_sec,overlap);
nwindows = height(windows);
edges = -1:.05:1;
results = zeros(nwindows,numel(edges)-1);
avg_corr = zeros(nwindows,1);
tic
for n = 1:nwindows
    disp(['Window #',num2str(n),'/',num2str(nwindows)]);
    thiswindow = windows(n,:);
    v.dataFilter = dft;
v.dataFilter.subjectGroup = 'naïve';
v.dataFilter.stims_allowed = 'all stimuli';
v.dataFilter.trial_sorting = 'chronological';
v.dataFilter.interval = [-22 -2];
[~,base_activity,labs] = ModeSelector(v).extract;
base_activity = cellfun(@(x) squeeze(mean(x,1,'omitmissing')),base_activity,'UniformOutput',false);
v.dataFilter.interval = thiswindow;
[~,odor_activity] = ModeSelector(v).extract;
odor_activity = cellfun(@(x) squeeze(mean(x,1,'omitmissing')),odor_activity,'UniformOutput',false);
labs = labs{1};
nsubjects = numel(base_activity);
ntrials = numel(labs);
y = zeros(ntrials,nsubjects);
for i = 1:nsubjects
    for j = 1:ntrials
        basepattern = base_activity{i}(:,j);
        odorpattern = odor_activity{i}(:,j);
        
        c = corrcoef(basepattern,odorpattern);
        y(j,i) = c(2);
    end
end
% t = 1:ntrials;
% mu = mean(y,2,'omitmissing');
% err = std(y,[],2,'omitmissing');
% figure;
% subplot(121)
% plotLineNShade(t,mu,err,'k',cfg);
% xticks(1:ntrials); xticklabels(labs)
% subplot(122)
% histogram(y(:),50,'FaceColor','k','EdgeAlpha',0)
% xlabel('r')
thisres = histcounts(y(:),edges);
results(n,:) = thisres(:);
avg_corr(n) = mean(y(:),'omitmissing');
end
toc

figure; imagesc(edges(1:end-1), windows(:,1),results)
hold on; plot(avg_corr,windows(:,1),'r-','LineWidth',cfg.lineWidth);
ylabel('Time from stim.onset (s)')
xlabel('r')
axis square; view([-90 90])
a = colorbar('Color',cfg.axcol);
a.Label.String = 'Count';
a.Label.FontSize= cfg.fontSize;
cfg.figSize = 'medium';
cfg.setFigure
cfg.saveFigure(gcf,'naive corr with prestimulus', saveType)
cfg = v.plotConfig;


%% naive corr with prestimulus in time
windows = defineTimeWindows(3,t_lim_sec,overlap);
nwindows = height(windows);
avg_corr = zeros(nwindows,nwindows,2);
tic
for n1 = 1:nwindows
for n2 = 1:nwindows
    disp(['Window #',num2str(n1),' vs #',num2str(n2),'/',num2str(nwindows)]);
    thiswindow1 = windows(n1,:);
    thiswindow2 = windows(n2,:);
    v.dataFilter = dft;
v.dataFilter.subjectGroup = 'trained';
v.dataFilter.stims_allowed = 'all stimuli';
v.dataFilter.trial_sorting = 'chronological';
v.dataFilter.interval = thiswindow2;
[~,base_activity,labs] = ModeSelector(v).extract;
base_activity = cellfun(@(x) squeeze(mean(x,1,'omitmissing')),base_activity,'UniformOutput',false);
v.dataFilter.interval = thiswindow1;
[~,odor_activity] = ModeSelector(v).extract;
odor_activity = cellfun(@(x) squeeze(mean(x,1,'omitmissing')),odor_activity,'UniformOutput',false);
labs = labs{1};
nsubjects = numel(base_activity);
ntrials = numel(labs);
y = zeros(ntrials,nsubjects);
for i = 1:nsubjects
    for j = 1:ntrials
        basepattern = base_activity{i}(:,j);
        odorpattern = odor_activity{i}(:,j);
        
        c = corrcoef(basepattern,odorpattern);
        y(j,i) = c(2);
    end
end
avg_corr(n1,n2,1) = mean(y(:),'omitmissing');
avg_corr(n1,n2,2) = std(y(:),'omitmissing');
end
end
toc

figure; imagesc(windows(:,1),windows(:,1),avg_corr(:,:,1))
xlabel('Time from stim.onset (s)')
ylabel('Time from stim.onset (s)')
axis square
a = colorbar('Color',cfg.axcol);
a.Label.String = 'r';
a.Label.FontSize= cfg.fontSize;
cfg.figSize = 'medium';
cfg.setFigure
cfg.saveFigure(gcf,'naive corr with prestimulus', saveType)
cfg = v.plotConfig;


%%
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'trained';
v.dataFilter.stims_allowed = 'all stimuli';
v.dataFilter.interval = [-5 35];
v.dataFilter.repetitions = [1:5];
[~,events,all_labs] = ModeSelector(v).extract;
% events = cellfun(@(x) movmean(x,3,1,'omitmissing'),events,'UniformOutput',false);
events = cellfun(@(x) permute(x,[3,2,1]),events,'UniformOutput',false);
events = separateTrials(events);
T = size(events{1},3);
figure; C = plotDistances(events,'full','correlation',1:T,cfg);
clim([0 1])
xlabel('Time from stim.onset (s)')
ylabel('Time from stim.onset (s)')

function out = separateTrials(data)
n = numel(data);
out = {};
for i = 1:n
    ntrials = size(data{i},1);
    for j = 1:ntrials
        out{end+1} = data{i}(j,:,:);
    end
end
out = out(:);
end



%% naive average response trace
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'naïve';
v.dataFilter.interval = [-5 40];
v.plotAvgResponseTrace;
v.dataFilter = dft;
cfg.setFigure
ylim([0 .2])
cfg.saveFigure(gcf,'naive average response trace', saveType)
cfg = v.plotConfig;

%% naive individual response traces (pooled imagesc)
Nshow = 1000; % n of cells to randomly subsample for visualization
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'naïve';
v.dataFilter.interval  = [.5 20]; % get intensity during odor window
iFR_odor = v.plotUnitActivityMetricHead('method','avg intensity');
N = cellfun(@height,iFR_odor); N = cumsum(N);
iFR_odor = cell2mat(iFR_odor);
v.dataFilter.interval = []; % get cellwise mean and std from all data
iFR_var = v.plotUnitActivityMetricHead('method','variance');
iFR_var = cell2mat(iFR_var);
iFR_std = sqrt(iFR_var); % variance to std
iFR_mu = v.plotUnitActivityMetricHead('method','avg intensity');
iFR_mu = cell2mat(iFR_mu);
iFR_odor = mean(iFR_odor,2,'omitmissing'); % avg over trials
iFR_std = mean(iFR_std,2,'omitmissing'); % avg over trials
iFR_mu = mean(iFR_mu,2,'omitmissing'); % avg over trials
resp = (iFR_odor-iFR_mu)./iFR_std; % responsiveness score
% plot responsiveness distributions
nsubjects = numel(N);
resp_all = cell(nsubjects,1);
N = [1;N];
figure; 
for i = 1:nsubjects
    resp_all{i} = resp(N(i):N(i+1),:,:);
    cdfplot(resp_all{i}); hold on
end
xlabel('Responsiveness (STD)')
ylabel('CDF'); title('')
axis square
cfg.setFigure
cfg.saveFigure(gcf,'naive responsiveness cdf', saveType)
cfg = v.plotConfig;
% plot individual responses
[resp, idx] = sort(resp,'descend');
v.dataFilter.interval = [-5 35]; % visualization windows
[~, events] = ModeSelector(v).extract;
events = cellfun(@(x) permute(x,[2 1 3]), events, 'UniformOutput', false);
events = cell2mat(events); % [N, T, trials]
[N, L, ntrials] = size(events);
t = linspace(v.dataFilter.interval(1),v.dataFilter.interval(2),L);
events = events(idx,:,:); % sort by responsiveness
idx = false(N,1); idx(randperm(N,Nshow)) = true; % random subsample
y = mean(events(idx,:,:),3,'omitmissing')-mean(events(idx,1:38,:),[2,3],'omitmissing');
figure; imagesc(t,1:Nshow,y)
clim([-0.2 .5])
xlabel('Time from stim. onset (s)')
ylabel('Cell #')
a = colorbar('Color',cfg.axcol);
a.Label.String = 'delta iFR (Hz)';
cfg.figSize = "medium"; cfg.aspRatioType = "tall";
cfg.setFigure
cfg.saveFigure(gcf,'naive subsample individual responses', saveType)
cfg = v.plotConfig;
figure; plot(resp(idx),'k');
ylim([min(ylim), 3])
xlabel('Cell #')
ylabel('Responsiveness (STD)')
title('Random subsample of cells')
cfg.setFigure
cfg.saveFigure(gcf,'naive subsample responsiveness', saveType)
cfg = v.plotConfig;
v.dataFilter = dft;

%% firing rate distributions
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'naïve';
v.dataFilter.traceType = 'pSpike';
v.dataFilter.interval  = []; % all data
[~,iFR_odor] = ModeSelector(v).extract;
iFR_odor = cellfun(@(x) permute(x,[2,1,3]),iFR_odor,'UniformOutput',false);
nsubjects = numel(iFR_odor);
edges = linspace(0,10,101);
Y = [];
for i = 1:nsubjects
    n = numel(iFR_odor{i});
    y = histcounts(iFR_odor{i},edges);
    Y = [Y; y(:)'./n];
end
figure;
hold on
plot(edges(1:end-1),Y,'k-','LineWidth',1)
xscale log; yscale log
xlabel('iFR (Hz)')
ylabel('Portion of frames')
axis square
cfg.setFigure
cfg.saveFigure(gcf,'naive firing rate distributions', saveType)
cfg = v.plotConfig;
v.dataFilter = dft;

%% 
close all