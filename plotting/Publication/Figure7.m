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
[group_data, ~] = loadAndPrepareGCMCData(v, folder_tag);

% --- All vs All odor manifolds ---
subdir = fullfiletol(gcmc_savepath, 'allvsall');
if ~isfolder(subdir); mkdir(subdir); end
f = GCMCResultsFilter('shuffle', false);
plotAndSaveBoxplots(group_data, cfg, f, metrics_to_plot, metric_labels, metric_yranges, subdir, saveType);

% --- Familiar odor manifolds ---
subdir = fullfiletol(gcmc_savepath, 'familiar');
if ~isfolder(subdir); mkdir(subdir); end
f = GCMCResultsFilter('shuffle', false, 'manifold_name_1', familiar_stims, 'manifold_name_2', familiar_stims);
plotAndSaveBoxplots(group_data, cfg, f, metrics_to_plot, metric_labels, metric_yranges, subdir, saveType);

% --- Novel odor manifolds ---
subdir = fullfiletol(gcmc_savepath, 'novel');
if ~isfolder(subdir); mkdir(subdir); end
f = GCMCResultsFilter('shuffle', false, 'manifold_name_1', novel_stims, 'manifold_name_2', novel_stims);
plotAndSaveBoxplots(group_data, cfg, f, metrics_to_plot, metric_labels, metric_yranges, subdir, saveType);

% ========== TRIAL MANIFOLDS (line plots by rep) ==========
folder_tag = 'trials';
gcmc_savepath = fullfiletol(savepath, 'gcmc', folder_tag);
if ~isfolder(gcmc_savepath); mkdir(gcmc_savepath); end

% Load data once
[group_data_trials, ~] = loadAndPrepareGCMCData(v, folder_tag);

% For rep-breakdown, we need to restructure group_data by (group, rep)
% But first, let's do line plots with the new plotting function

% --- Same rep trial manifolds (all stimuli) ---
subdir = fullfiletol(gcmc_savepath, 'allstims_samerep');
if ~isfolder(subdir); mkdir(subdir); end
f = GCMCResultsFilter('shuffle', false);
f.pairFilters.manifold_name = 'same';  % same stimulus only
f.pairFilters.manifold_rep = 'same';   % same repetition only
plotAndSaveRepLines(group_data_trials, cfg, f, metrics_to_plot, metric_labels, subdir, saveType);

% --- Same rep trial manifolds (familiar stimuli only) ---
subdir = fullfiletol(gcmc_savepath, 'familiar_samerep');
if ~isfolder(subdir); mkdir(subdir); end
f = GCMCResultsFilter('shuffle', false, 'manifold_name_1', familiar_stims, 'manifold_name_2', familiar_stims);
f.pairFilters.manifold_name = 'same';
f.pairFilters.manifold_rep = 'same';
plotAndSaveRepLines(group_data_trials, cfg, f, metrics_to_plot, metric_labels, subdir, saveType);

% --- Same rep trial manifolds (novel stimuli only) ---
subdir = fullfiletol(gcmc_savepath, 'novel_samerep');
if ~isfolder(subdir); mkdir(subdir); end
f = GCMCResultsFilter('shuffle', false, 'manifold_name_1', novel_stims, 'manifold_name_2', novel_stims);
f.pairFilters.manifold_name = 'same';
f.pairFilters.manifold_rep = 'same';
plotAndSaveRepLines(group_data_trials, cfg, f, metrics_to_plot, metric_labels, subdir, saveType);



function [group_data, all_results] = loadAndPrepareGCMCData(v, folder_tag)
    % Load GCMC results and prepare group data
    indir = fullfiletol('manifold_data', folder_tag);
    all_results = GCMC_Analysis(v).extractResultsFromMultipleSubjects(fullfiletol(indir, 'odorexp004_IC1_130625'));
    group_data_raw = GCMC_Analysis(v).clusterByGroup(all_results);
    group_data = mergeTrainedGroups(group_data_raw);
end

function new_group_data = mergeTrainedGroups(group_data)
    % Merge trained subgroups into single 'trained' group
    % Expected order: naïve, trained1, trained2, trained1(...), uncoupled
    new_group_data = struct('group_name', {}, 'data', {});

    % Find naïve and uncoupled
    naive_idx = find(strcmp({group_data.group_name}, 'naïve'));
    uncoupled_idx = find(strcmp({group_data.group_name}, 'uncoupled'));
    trained_idx = find(contains({group_data.group_name}, 'trained'));

    idx = 0;
    if ~isempty(naive_idx)
        idx = idx + 1;
        new_group_data(idx) = group_data(naive_idx);
    end
    if ~isempty(trained_idx)
        idx = idx + 1;
        new_group_data(idx).group_name = 'trained';
        new_group_data(idx).data = group_data(trained_idx(1)).data;
        for i = 2:numel(trained_idx)
            new_group_data(idx).data = [new_group_data(idx).data; group_data(trained_idx(i)).data];
        end
    end
    if ~isempty(uncoupled_idx)
        idx = idx + 1;
        new_group_data(idx) = group_data(uncoupled_idx);
    end
end

function plotAndSaveBoxplots(group_data, cfg, filter, metrics, labels, yranges, savedir, saveType)
    % Plot and save boxplots for multiple metrics
    % plotBoxplotsByGroup returns array of figure handles (one per metric it plots)

    hf = GCMC_Plotting.plotBoxplotsByGroup(group_data, cfg, filter);

    % The function plots 7 metrics by default (hardcoded in plotBoxplotsByGroup)
    % We only want to save our specified metrics
    nFigs = min(numel(hf), numel(metrics));

    for i = 1:nFigs
        if i <= numel(labels) && i <= numel(yranges)
            figure(hf(i));
            ylim(yranges{i});
            title('');
            ylabel(labels{i});
            cfg.figSize = 'tiny';
            cfg.aspRatioType = 'tall';
            cfg.lineWidth = 0.5;
            cfg.setLines = true;
            cfg.setFigure;
            cfg.savePath = savedir;
            cfg.saveFigure(gcf, [labels{i}, ' boxplot'], saveType);
        end
    end
    close all;
end

function plotAndSaveRepLines(group_data, cfg, filter, metrics, labels, savedir, saveType)
    % Plot and save line plots across repetitions for multiple metrics
    for i = 1:numel(metrics)
        [hf, stats] = GCMC_Plotting.plotMetricsByRepLine(group_data, cfg, filter, metrics{i});

        % Display stats summary
        disp(['=== Stats for ', labels{i}, ' ===']);
        for j = 1:numel(stats.friedman)
            disp(sprintf('  Friedman %s: p=%.4g', stats.friedman(j).group, stats.friedman(j).p));
        end
        sig_between = find([stats.between_group.p_fdr] < 0.05);
        for j = sig_between
            disp(sprintf('  Between-group rep%d %s vs %s: p_fdr=%.4g', ...
                stats.between_group(j).rep, stats.between_group(j).group1, ...
                stats.between_group(j).group2, stats.between_group(j).p_fdr));
        end

        ylabel(labels{i});
        title('');
        cfg.figSize = 'small';
        cfg.aspRatioType = 'wide';
        cfg.setFigure;
        cfg.savePath = savedir;
        cfg.saveFigure(gcf, [labels{i}, ' by rep'], saveType);
        close(hf);
    end
end


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