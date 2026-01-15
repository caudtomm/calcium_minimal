
dbstop if error

s = true; % save figures to files?
savepath = 'bin2';
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
cfg = PlotConfig;
cfg.savePath = savepath;

v = ExperimentViewer(experiment);
v.plotConfig = cfg;
cfg = v.plotConfig;

v.dataFilter.traceType = 'pSpike';
v.dataFilter.interval = [1,20];
v.dataFilter.trial_sorting = 'stim_id';
dft = v.dataFilter;
 

%% FIGURE 1

%% naive average response trace
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'naïve';
v.dataFilter.interval = [-5 40];
v.plotAvgResponseTrace;
v.dataFilter = dft;
cfg.setFigure
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