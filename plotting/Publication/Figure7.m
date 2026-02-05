dbstop if error

s = true; % save figures to files?
savepath = 'bin7';
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
cfg = PlotConfig('colormapName','lapaz','favouriteColors',[84,85,73,86:99]); % (test1, test2, ctrl)
cfg.custom.crange = [.3 .7];
cfg.savePath = savepath;

v = ExperimentViewer(experiment);
v.plotConfig = cfg;
cfg = v.plotConfig;

v.dataFilter.traceType = 'pSpike';
v.dataFilter.interval = [1,20];
v.dataFilter.trial_sorting = 'stim_id';
dft = v.dataFilter;
 

%% FIGURE 7


%% plot NMF lines
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'trained';
v.dataFilter.stims_allowed = 'all stimuli';
v.dataFilter.interval = [.5, 20];
v.dataFilter.repetitions = 1:5;
v.dataFilter.mode_name = 'nmf';

[~,events,labs] = ModeSelector(v).extract;
proj = computeLDE(events,labs,'pooldata',true,'nans2zeros',true, 'method','returnonly');
figure; out = plotLDE(proj.embedding{1}.reduction(:,1:3,:),'lines',proj.labs{1},cfg, 'ldetype','nmf'); % plot
% plotting the first 3 NMF components


%% plot PC variance explained for all groups
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'all';
v.dataFilter.stims_allowed = 'all stimuli';
v.dataFilter.interval = [.5, 20];
v.dataFilter.repetitions = 1:5;

[~,events,labs] = ModeSelector(v).extract;
proj = computeLDE(events,labs,'pooldata',true,'nans2zeros',true, 'method','pca');

emb = proj.embedding{1};
dt = proj.inputData{1};
out = characterizePCspace(dt.traces,emb.coeff,emb.explained,cfg,50);
xscale log; yscale log;
xlim([min(xlim) out.maxPC+100])
cfg.setFigure;
cfg.saveFigure(gcf,'allgroups PC variance explained', saveType)


%% plot PCA lines
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'naïve';
v.dataFilter.stims_allowed = 'all stimuli';
v.dataFilter.interval = [.5, 20];
v.dataFilter.repetitions = 1:5;

[~,events,labs] = ModeSelector(v).extract;
proj = computeLDE(events,labs,'pooldata',true,'nans2zeros',true, 'method','pca');
figure; out = plotLDE(proj.embedding{1}.reduction(:,1:2,:),'lines',proj.labs{1},cfg, 'ldetype',proj.name); % plot
xlabel('PC 1'); ylabel('PC 2');
cfg.setLines = false;
cfg.figSize = 'medium';
cfg.setFigure;
cfg.saveFigure(gcf,'naive PCA 2d', saveType)


%% plot UMAP lines
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'uncoupled';
v.dataFilter.stims_allowed = 'all stimuli';
v.dataFilter.interval = [.5, 20];
v.dataFilter.repetitions = 1:5;

[~,events,labs] = ModeSelector(v).extract;
proj = computeLDE(events,labs,'pooldata',true,'nans2zeros',true, 'method','umap');
figure; out = plotLDE(proj.embedding{1}.reduction,...
    'lines',proj.labs{1},cfg, 'ldetype',proj.name); % plot
xlabel('UMAP 1'); ylabel('UMAP 2');
cfg.setLines = false;
cfg.figSize = 'medium';
cfg.setFigure;
cfg.saveFigure(gcf,'uncoupled UMAP 2d', saveType)


%%

grouptag = 'naive';

v.dataFilter = dft; % recover
v.dataFilter.subjectGroup = 'naïve';
R = driftMetricsFigure(v);

figure(R.figures.angle_box)
axis normal; ylim([0 1]); xlim([0 max(xlim)])
xtickangle(90)
cfg.figSize = 'tiny';
cfg.aspRatioType = 'tall';
cfg.setLines = false;
cfg.setFigure;
cfg.saveFigure(gcf,[grouptag, ' interrep angles boxplot'], saveType)

figure(R.figures.anglecorr_box)
axis normal; ylim([0 1]); xlim([0 max(xlim)])
xtickangle(90)
cfg.figSize = 'tiny';
cfg.aspRatioType = 'tall';
cfg.setLines = false;
cfg.setFigure;
cfg.saveFigure(gcf,[grouptag, ' interrep angle cosine boxplot'], saveType)

figure(R.figures.corr_box)
axis normal; ylim([-1 1]); xlim([0 max(xlim)])
xtickangle(90)
cfg.figSize = 'tiny';
cfg.aspRatioType = 'tall';
cfg.setLines = false;
cfg.setFigure;
cfg.saveFigure(gcf,[grouptag, ' interrep vector corr boxplot'], saveType)


%% GCMC
close all

% Metrics to plot
metrics_to_plot = {'Fun_capacity', 'Fun_dimension', 'Fun_radius', ...
                   'Fun_center_alignment', 'Fun_axis_alignment'};
metric_labels = {'Capacity', 'Dimension', 'Radius', 'Center alignment', 'Axis alignment'};
metric_yranges = {[0 .15], [0 30], [0 2], [0 1], [0 .5]};

% Stimulus sets
familiar_stims = {'Arg', 'Ala', 'His'};
novel_stims = {'Trp', 'Ser', 'Leu'};

%% ========== ODOR MANIFOLDS (boxplots by group) ==========
folder_tag = 'odors';
gcmc_savepath = fullfiletol(savepath, 'gcmc', folder_tag);
if ~isfolder(gcmc_savepath); mkdir(gcmc_savepath); end

% Load data once
exp_name = 'odorexp004_IC1_130625';
[group_data, ~] = GCMC_Plotting.loadAndPrepareData(v, folder_tag, exp_name);

% --- All vs All odor manifolds ---
subdir = fullfiletol(gcmc_savepath, 'allvsall');
if ~isfolder(subdir); mkdir(subdir); end
f = GCMCResultsFilter('shuffle', false);
GCMC_Plotting.plotAndSaveBoxplots(group_data, cfg, f, metrics_to_plot, metric_labels, metric_yranges, subdir, saveType);

% --- Familiar odor manifolds ---
subdir = fullfiletol(gcmc_savepath, 'familiar');
if ~isfolder(subdir); mkdir(subdir); end
f = GCMCResultsFilter('shuffle', false, 'manifold_name_1', familiar_stims, 'manifold_name_2', familiar_stims);
GCMC_Plotting.plotAndSaveBoxplots(group_data, cfg, f, metrics_to_plot, metric_labels, metric_yranges, subdir, saveType);

% --- Novel odor manifolds ---
subdir = fullfiletol(gcmc_savepath, 'novel');
if ~isfolder(subdir); mkdir(subdir); end
f = GCMCResultsFilter('shuffle', false, 'manifold_name_1', novel_stims, 'manifold_name_2', novel_stims);
GCMC_Plotting.plotAndSaveBoxplots(group_data, cfg, f, metrics_to_plot, metric_labels, metric_yranges, subdir, saveType);

%% ========== TRIAL MANIFOLDS (line plots by rep) ==========
folder_tag = 'trials_';
gcmc_savepath = fullfiletol(savepath, 'gcmc', folder_tag);
if ~isfolder(gcmc_savepath); mkdir(gcmc_savepath); end

% Load data once
[group_data_trials, ~] = GCMC_Plotting.loadAndPrepareData(v, folder_tag, exp_name);

% --- Same rep trial manifolds (all stimuli) ---
subdir = fullfiletol(gcmc_savepath, 'allstims_samerep');
if ~isfolder(subdir); mkdir(subdir); end
f = GCMCResultsFilter('shuffle', false);
f.pairFilters.manifold_name = 'different';  % different stimulus only
f.pairFilters.manifold_rep = 'same';   % same repetition only
GCMC_Plotting.plotAndSaveRepLines(group_data_trials, cfg, f, metrics_to_plot, metric_labels, subdir, saveType);

% --- Same rep trial manifolds (familiar stimuli only) ---
subdir = fullfiletol(gcmc_savepath, 'familiar_samerep');
if ~isfolder(subdir); mkdir(subdir); end
f = GCMCResultsFilter('shuffle', false, 'manifold_name_1', familiar_stims, 'manifold_name_2', familiar_stims);
f.pairFilters.manifold_name = 'different';
f.pairFilters.manifold_rep = 'same';
GCMC_Plotting.plotAndSaveRepLines(group_data_trials, cfg, f, metrics_to_plot, metric_labels, subdir, saveType);

% --- Same rep trial manifolds (novel stimuli only) ---
subdir = fullfiletol(gcmc_savepath, 'novel_samerep');
if ~isfolder(subdir); mkdir(subdir); end
f = GCMCResultsFilter('shuffle', false, 'manifold_name_1', novel_stims, 'manifold_name_2', novel_stims);
f.pairFilters.manifold_name = 'different';
f.pairFilters.manifold_rep = 'same';
GCMC_Plotting.plotAndSaveRepLines(group_data_trials, cfg, f, metrics_to_plot, metric_labels, subdir, saveType);


%% template matching lines
classifier = 'qda';
v.dataFilter = dft;
grouptag = 'naive';
stimtag = 'allstims';
v.dataFilter.subjectGroup = 'naïve';
v.dataFilter.stims_allowed = 'all stimuli';
v.dataFilter.interval = [.5, 20];
v.dataFilter.repetitions = 1:5;
method = 'correlation';
focus_stims = 'all stimuli';

hf = figure;
c = v.plotDiscriminationHead(...
    'plotType', 'performance_lines', ...
    'method',method, ...
    'focus_stims', focus_stims, ...
    'classifier',classifier,...
    'trainblockmode','2blocks',...
    'separatetestset',true,...
    'zscore',false);
for i = 2:width(c)
    p = signrank(c(:,1),c(:,i));
    disp(['Paired Wilcoxon signed-rank test - reps 1 vs ', num2str(i),': ',num2str(p)])
end
xlabel('Template trial #')
ylabel('Performance')
ylim([0 1])
title(classifier)
cfg.figSize = "tiny";
cfg.aspRatioType = "tall";
cfg.setFigure;
cfg.saveFigure(gcf,['template ',grouptag, ' ',stimtag,' lines'], saveType)


%% template matching mats

v.dataFilter = dft;
grouptag = 'uncoupled';
stimtag = 'allstims';
v.dataFilter.subjectGroup = 'uncoupled';
v.dataFilter.stims_allowed = 'all stimuli';
v.dataFilter.interval = [.5, 20];
v.dataFilter.repetitions = 1:5;
method = 'correlation';
focus_stims = 'all stimuli';
do_zscore = false;

hf = figure;
out = v.plotDiscriminationHead(...
    'plotType', 'performance_mat', ...
    'method',method, ...
    'focus_stims', focus_stims, ...
    'zscore',do_zscore);
clim(cfg.custom.crange)
xticks([]); yticks([])
cfg.figSize = "small";
cfg.aspRatioType = "square";
cfg.setFigure;
cfg.saveFigure(gcf,['template ',grouptag, ' ',stimtag,' mat'], saveType)