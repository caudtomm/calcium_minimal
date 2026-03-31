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
cfg.custom.crange = [.1 .7];
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
[out, ax1, ax2] = characterizePCspace(dt.traces,emb.coeff,emb.explained,cfg,50);
axes(ax1)
xscale log; yscale log;
xlim([min(xlim) out.maxPC+100])
axes(ax2)
xscale log; yscale log;
% xlim([min(xlim) out.maxPC+100])
axes(ax1)
ylim([0.01 max(ylim)*2])
cfg.setFigure;
cfg.saveFigure(gcf,'allgroups PC variance explained', saveType)


%% plot PCA lines
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'all';
v.dataFilter.stims_allowed = 'all stimuli';
v.dataFilter.interval = [.5, 20];
v.dataFilter.repetitions = 1:5;

[~,events,labs] = ModeSelector(v).extract;
proj = computeLDE(events,labs,'pooldata',true,'nans2zeros',true, 'method','pca');
figure; out = plotLDE(proj.embedding{1}.reduction(:,1:2,:),'lines',proj.labs{1},cfg, 'ldetype',proj.name); % plot
xlabel('PC 1'); ylabel('PC 2');
legend off
cfg.setLines = false;
cfg.figSize = 'small';
cfg.setFigure;
cfg.saveFigure(gcf,'allgroups PCA 2d raw', saveType)


%% plot UMAP lines
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'all';
v.dataFilter.stims_allowed = 'all stimuli';
v.dataFilter.interval = [1, 20];
v.dataFilter.repetitions = 1:5;

% v.dataFilter.mode_name = 'dpca';
% v.dataFilter.mode_OI = 'stimulus';
% v.dataFilter.mode_method = 'isolate';
% v.dataFilter.mode_file = 'dpca_trained.mat';

[~,events,labs] = ModeSelector(v).extract;

% [~,odor_events,all_labs] = ModeSelector(v).extract;
% L = height(odor_events{1});
% v.dataFilter.interval = [-22 -2];
% [~,base_events] = ModeSelector(v).extract;
% base_events = cellfun(@(x) repmat(mean(x,1,'omitmissing'),L,1,1),base_events,'UniformOutput',false);
% events = cellfun(@(x,y) y-x,base_events,odor_events,'UniformOutput',false);

proj = computeLDE(events,labs,'pooldata',true,'nans2zeros',true, ...
    'method','umap','n_components',2, 'metric','euclidean');
figure; out = plotLDE(proj.embedding{1}.reduction,...
    'lines',proj.labs{1},cfg, 'ldetype',proj.name); % plot
xlabel('UMAP 1'); ylabel('UMAP 2');
legend off
cfg.setLines = false;
cfg.figSize = 'small';
cfg.setFigure;
cfg.saveFigure(gcf,'allgroups UMAP 2d raw', saveType)


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
exp_name = 'odorexp004_IC1_130625';

% Metrics to plot
metrics_to_plot = {'Fun_capacity', 'Fun_dimension', 'Fun_radius', ...
                   'Fun_center_alignment', 'Fun_axis_alignment'};
metric_labels = {'Capacity', 'Dimension', 'Radius', 'Center alignment', 'Axis alignment'};
metric_yranges = {[0 .15], [0 30], [0 2], [0 1], [0 .5]}; % for odor manifolds
%metric_yranges = {[.3 .6], [3 4.5], [0.8 1.1], [.4 .7], [.1 .2]}; % for sliding windows

% Stimulus sets
familiar_stims = {'Arg', 'Ala', 'His'};
novel_stims = {'Trp', 'Ser', 'Leu'};

%% ========== ODOR MANIFOLDS (boxplots by group) ==========
folder_tag = 'repetitions';
gcmc_savepath = fullfiletol(savepath, 'gcmc', folder_tag);
if ~isfolder(gcmc_savepath); mkdir(gcmc_savepath); end

% Load data once
[group_data, ~] = GCMC_Plotting.loadAndPrepareData(v, folder_tag, exp_name);
group_data = group_data(1:2)

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
folder_tag = 'trials';
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
stats = GCMC_Plotting.plotAndSaveRepLines(group_data_trials, cfg, f, metrics_to_plot, metric_labels, subdir, saveType,'average_by_subject',true);

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


%% ========== SLIDING WINDOW TRIAL MANIFOLDS ==========
folder_tag = 'trial_slide_windows';
gcmc_savepath = fullfiletol(savepath, 'gcmc', folder_tag);
if ~isfolder(gcmc_savepath); mkdir(gcmc_savepath); end

% Load all time windows
sw = GCMC_Plotting.loadSlidingWindowData(v, folder_tag, exp_name);

% Color limits per metric: [nMetrics x 2] - same as metric_yranges but for heatmaps
sw_clim = cell2mat(metric_yranges');  % convert {[0 .15], [0 30], ...} to [nM x 2]

% --- All stimuli, same-rep, different-stim (absolute) ---
subdir = fullfiletol(gcmc_savepath, 'allstims_samerep');
if ~isfolder(subdir); mkdir(subdir); end
f = GCMCResultsFilter('shuffle', false);
f.pairFilters.manifold_name = 'different';
f.pairFilters.manifold_rep = 'same';
[hf_sw, pf_sw] = GCMC_Plotting.plotSlidingWindowMetrics(sw, cfg, f, metrics_to_plot, ...
    'clim', sw_clim,'relative',false);

% Display Friedman p-values
for i_m = 1:numel(metrics_to_plot)
    disp(['=== Sliding window Friedman: ', metric_labels{i_m}, ' ===']);
    disp(array2table(pf_sw(:,:,i_m), ...
        'VariableNames', sw.windows, 'RowNames', sw.groups));
end

% Save figures
for i = 1:numel(hf_sw)
    if isgraphics(hf_sw(i))
        figure(hf_sw(i));
        cfg.savePath = subdir;
        cfg.saveFigure(hf_sw(i), get(get(gca,'Title'),'String'), saveType);
    end
end
close all; clear hf_sw

% Plot p-values as lines and save
for i_g = 1:numel(sw.groups)
    for i_m = 1:numel(metrics_to_plot)
        thisdt = pf_sw(i_g,:,i_m);
        hf = figure; plot(thisdt,'k-');
        hold on; yline(0.05,'r-')
        yscale log
        xticks(1:numel(sw.windows)); xticklabels(sw.t_centers)
        ylabel('Friedman p-value')
        ttlstr = [sw.groups{i_g},' - ' metric_labels{i_m}];
        title(ttlstr)
        cfg.savePath = subdir;
        cfg.setFigure
        cfg.saveFigure(gcf, ['pval - ',ttlstr], saveType);
        close
    end
end

% --- All stimuli, same-rep, different-stim (RELATIVE to pre-stimulus) ---
subdir = fullfiletol(gcmc_savepath, 'allstims_samerep_rel');
if ~isfolder(subdir); mkdir(subdir); end
[hf_sw, ~] = GCMC_Plotting.plotSlidingWindowMetrics(sw, cfg, f, metrics_to_plot, ...
    'relative', true);
for i = 1:numel(hf_sw)
    if isgraphics(hf_sw(i))
        figure(hf_sw(i));
        cfg.savePath = subdir;
        cfg.saveFigure(hf_sw(i), get(get(gca,'Title'),'String'), saveType);
    end
end
close all;

% --- Familiar stimuli, same-rep ---
subdir = fullfiletol(gcmc_savepath, 'familiar_samerep');
if ~isfolder(subdir); mkdir(subdir); end
f = GCMCResultsFilter('shuffle', false, 'manifold_name_1', familiar_stims, 'manifold_name_2', familiar_stims);
f.pairFilters.manifold_name = 'different';
f.pairFilters.manifold_rep = 'same';
[hf_sw, ~] = GCMC_Plotting.plotSlidingWindowMetrics(sw, cfg, f, metrics_to_plot, ...
    'clim', sw_clim);
for i = 1:numel(hf_sw)
    if isgraphics(hf_sw(i))
        figure(hf_sw(i));
        cfg.savePath = subdir;
        cfg.saveFigure(hf_sw(i), get(get(gca,'Title'),'String'), saveType);
    end
end
close all;

% --- Novel stimuli, same-rep ---
subdir = fullfiletol(gcmc_savepath, 'novel_samerep');
if ~isfolder(subdir); mkdir(subdir); end
f = GCMCResultsFilter('shuffle', false, 'manifold_name_1', novel_stims, 'manifold_name_2', novel_stims);
f.pairFilters.manifold_name = 'different';
f.pairFilters.manifold_rep = 'same';
[hf_sw, ~] = GCMC_Plotting.plotSlidingWindowMetrics(sw, cfg, f, metrics_to_plot, ...
    'clim', sw_clim);
for i = 1:numel(hf_sw)
    if isgraphics(hf_sw(i))
        figure(hf_sw(i));
        cfg.savePath = subdir;
        cfg.saveFigure(hf_sw(i), get(get(gca,'Title'),'String'), saveType);
    end
end
close all;


%% decoding lines
classifier = 'template_match';
v.dataFilter = dft;
grouptag = 'naive';
stimtag = 'all';
v.dataFilter.subjectGroup = 'naïve';
v.dataFilter.stims_allowed = 'all stimuli'; %{'Trp','Ser','Leu'};
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
    'trainblockmode','single',...
    'separatetestset',true,...
    'zscore',false);
for i = 2:width(c)
    p = signrank(c(:,1),c(:,i));
    disp(['Paired Wilcoxon signed-rank test - reps 1 vs ', num2str(i),': ',num2str(p)])
end
fr_trained = statsUtils.friedman(c);
fprintf('Friedman test: chi2=%.2f, df=%d, p=%.4g\n', ...
    fr_trained.chi2, fr_trained.df, fr_trained.p);
xlabel('Template trial #')
ylabel('Performance')
ylim([0 1])
title(classifier)
legend off
cfg.figSize = "tiny";
cfg.aspRatioType = "tall";
cfg.setFigure;
cfg.saveFigure(gcf,[classifier,' ',grouptag, ' ',stimtag,' lines'], saveType)


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

%%

v.plotConfig.theme = 'many colors';
v.dataFilter.subjectGroup = 'all';
[R, ira] = driftMetricsFigure(v);
close all
angdiff_sham = getAngleDiff(ira);

function y = getAngleDiff(ira)
    [nstims, nangles, nsubj] = size(ira);
    x = nchoosek(1:nstims,2);
    npairs = height(x);
    y = [];
    for i = 1:npairs
        dt = ira(x(i,:),:,:); % [odor1;odor2 x angles x subj]
        dtd = diff(dt,[],1);
        y = [y; dtd(:)];
    end
end