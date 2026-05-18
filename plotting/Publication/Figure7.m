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
v.dataFilter.subjectGroup = 'trained';
v.dataFilter.stims_allowed = 'all stimuli';
v.dataFilter.trial_sorting = 'chronological';
v.dataFilter.interval = [.5 20];
v.dataFilter.repetitions = 1:5;
cfg.useNatureColors = true;

% v.dataFilter.mode_name = 'dpca';
% v.dataFilter.mode_OI = 'stimulus';
% v.dataFilter.mode_method = 'isolate';
% v.dataFilter.mode_file = 'dpca_naive_BS.mat';

% [~,events,labs] = ModeSelector(v).extract;
% events = filterEventsByScore(events,v,[.1 .9],'linear');

[~,odor_events,all_labs] = ModeSelector(v).extract;
L = height(odor_events{1});
v.dataFilter.interval = [-22 -2];
[~,base_events] = ModeSelector(v).extract;
base_events = cellfun(@(x) repmat(mean(x,1,'omitmissing'),L,1,1),base_events,'UniformOutput',false);
events = cellfun(@(x,y) y-x,base_events,odor_events,'UniformOutput',false);

proj = computeLDE(events,labs,'pooldata',true,'nans2zeros',true, 'method','pca');
figure; out = plotLDE(proj.embedding{1}.reduction(:,1:2,:),'trialNum',proj.labs{1},cfg, 'ldetype',proj.name); % plot
xlabel('PC 1'); ylabel('PC 2');
legend off
cfg.setLines = false;
cfg.figSize = 'small';
cfg.setFigure;
cfg.saveFigure(gcf,'trained BS PCA trialNum 2d raw', saveType)


%% plot UMAP lines
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'naïve';
v.dataFilter.stims_allowed = 'all stimuli';
v.dataFilter.interval = [1, 20];
v.dataFilter.repetitions = 1:5;
cfg.useNatureColors = true;

% v.dataFilter.mode_name = 'dpca';
% v.dataFilter.mode_OI = 'stimulus';
% v.dataFilter.mode_method = 'isolate';
% v.dataFilter.mode_file = 'dpca_naive_BS.mat';

[~,events,labs] = ModeSelector(v).extract;
events = filterEventsByScore(events,v,[.1 .9],'linear');

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
cfg.saveFigure(gcf,'naive low-drifters UMAP 2d raw', saveType)


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

%% Cell-level scatter: odor responsiveness vs baseline FR / corr. contribution / identity dPC score
groups_sc     = {'naïve', 'trained', 'uncoupled'};
dpca_files_sc = {'dpca_naive.mat', 'dpca_trained.mat', 'dpca_uncoupled.mat'};
fav_idx       = [84, 85, 73];   % lapaz indices for naive / trained / uncoupled

resp_sc = []; baseline_sc = []; cm_sc = []; mx_sc = []; dpc_sc = [];
gidx_sc = [];

for g = 1:numel(groups_sc)
    v.dataFilter = dft;
    v.dataFilter.subjectGroup  = groups_sc{g};
    v.dataFilter.stims_allowed = 'all stimuli';
    v.dataFilter.repetitions   = 1:5;

    % odor responsiveness: (mean_odor_FR - session_mean) / session_std
    v.dataFilter.interval = [1, 20];
    ifr_odor = v.plotUnitActivityMetricHead('method', 'avg intensity');
    ifr_odor = cell2mat(cellfun(@(x) mean(x,2,'omitmissing'), ifr_odor, 'UniformOutput',false));
    v.dataFilter.interval = [];
    ifr_var  = v.plotUnitActivityMetricHead('method', 'variance');
    ifr_mu   = v.plotUnitActivityMetricHead('method', 'avg intensity');
    ifr_std  = cell2mat(cellfun(@(x) sqrt(mean(x,2,'omitmissing')), ifr_var, 'UniformOutput',false));
    ifr_mu   = cell2mat(cellfun(@(x) mean(x,2,'omitmissing'), ifr_mu, 'UniformOutput',false));
    resp_g   = (ifr_odor - ifr_mu) ./ ifr_std;

    % mean baseline FR (pre-stimulus window)
    v.dataFilter.interval = [-22, -2];
    ifr_base = v.plotUnitActivityMetricHead('method', 'avg intensity');
    base_g   = cell2mat(cellfun(@(x) mean(x,2,'omitmissing'), ifr_base, 'UniformOutput',false));

    % per-cell contribution to inter-stimulus correlations: c_i = x_i*y_i / (||x||*||y||)
    v.dataFilter.interval    = [1, 20];
    v.dataFilter.repetitions = 1:5;
    tc_g = v.plotUnitActivityMetricHead('method', 'tuning curves');  % {sj}: [N x nS x nR]
    cmean_g = []; cmax_g = [];
    for sj = 1:numel(tc_g)
        tc = tc_g{sj};  [N, nS, nR] = size(tc);
        if nS < 2
            cmean_g = [cmean_g; nan(N,1)];
            cmax_g  = [cmax_g;  nan(N,1)];
            continue
        end
        pairs = nchoosek(1:nS, 2);
        C     = nan(N, size(pairs,1) * nR);
        col   = 0;
        for r = 1:nR
            for p = 1:size(pairs,1)
                col = col + 1;
                x = tc(:, pairs(p,1), r);
                y = tc(:, pairs(p,2), r);
                % d = sqrt(nansum(x.^2) * nansum(y.^2));
                % if d > eps; C(:,col) = x .* y / d; end
                xmean = mean(x,'omitmissing');
                ymean = mean(y,'omitmissing');
                C(:,col) = (x-xmean).*(y-ymean);
            end
        end
        cmean_g = [cmean_g; mean(C, 2, 'omitmissing')];
        cmax_g  = [cmax_g;  max(C, [], 2)];
    end

    % mean abs identity dPC loading
    v.dataFilter.mode_name   = 'dpca';
    v.dataFilter.mode_OI     = 'stimulus';
    v.dataFilter.mode_method = 'mode_values';
    v.dataFilter.mode_file   = dpca_files_sc{g};
    v.dataFilter.interval    = [1, 20];
    m_d   = ModeSelector(v).extract;
    dpc_g = nan(numel(resp_g), 1);
    ptr   = 0;
    for i = 1:numel(m_d.coeffs)
        if isempty(m_d.coeffs{i}); continue; end
        ids   = m_d.parseModeOI(i);
        W_id  = m_d.coeffs{i}(:, ids);   % [N x n_identity_modes]
        n     = size(W_id, 1);
        if ptr + n > numel(dpc_g); break; end
        dpc_g(ptr + (1:n)) = mean(abs(W_id), 2);
        ptr = ptr + n;
    end

    resp_sc     = [resp_sc;     resp_g];
    baseline_sc = [baseline_sc; base_g];
    cm_sc       = [cm_sc;       cmean_g];
    mx_sc       = [mx_sc;       cmax_g];
    dpc_sc      = [dpc_sc;      dpc_g];
    gidx_sc     = [gidx_sc;     repmat(g, numel(resp_g), 1)];
end

gcols_sc = cfg.c(fav_idx, :);   % [3 x 3] group colors

ydata_sc   = {baseline_sc,       cm_sc,                      mx_sc,                     dpc_sc};
col_titles = {'Mean baseline FR', 'Mean corr. contribution', 'Max corr. contribution',  'Mean |identity dPC|'};

hf_sc = figure;
tl = tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
xlabel(tl, 'Odor responsiveness (z)');

ax_sc = gobjects(3, 4);
for g = 1:3
    for k = 1:4
        ax_sc(g,k) = nexttile((g-1)*4 + k);
        mask = gidx_sc == g;
        scatter(ax_sc(g,k), resp_sc(mask), ydata_sc{k}(mask), 10, ...
            'filled', 'MarkerFaceAlpha', 0.4, 'MarkerFaceColor', gcols_sc(g,:));
        if g == 1; title(ax_sc(g,k), col_titles{k}); end
        if k == 1; ylabel(ax_sc(g,k), groups_sc{g}); end
        axis(ax_sc(g,k), 'square');  box(ax_sc(g,k), 'off');
        set(ax_sc(g,k), 'color', cfg.bgcol, 'XColor', cfg.axcol, 'YColor', cfg.axcol);
    end
end

% equalize y-limits per column across groups
for k = 1:4
    yl = cell2mat(arrayfun(@(g) ylim(ax_sc(g,k)), (1:3)', 'UniformOutput', false));
    shared_ylim = [min(yl(:,1)), max(yl(:,2))];
    for g = 1:3; ylim(ax_sc(g,k), shared_ylim); end
end

set(hf_sc, 'color', cfg.bgcol);
cfg.figSize = 'large';  cfg.aspRatioType = 'tall';  cfg.setLines = false;
cfg.setFigure;
cfg.saveFigure(hf_sc, 'odor resp vs cell metrics scatter', saveType)

%% Scatter correlation statistics: r and R² per subplot
stat_rows = cell(numel(groups_sc) * numel(col_titles), 6);
row = 0;
for g = 1:numel(groups_sc)
    for k = 1:numel(col_titles)
        row = row + 1;
        mask = gidx_sc == g;
        x = resp_sc(mask);
        y = ydata_sc{k}(mask);
        ok = isfinite(x) & isfinite(y);
        n = sum(ok);
        if n >= 3
            [r, p] = corr(x(ok), y(ok));
        else
            r = NaN;  p = NaN;
        end
        stat_rows(row,:) = {groups_sc{g}, col_titles{k}, n, r, r^2, p};
    end
end
corr_tab = cell2table(stat_rows, ...
    'VariableNames', {'Group', 'Metric', 'n', 'r', 'R2', 'p'});
disp(corr_tab)


%% PC1 trajectory vs trial number
group_pc1 = 'naïve';

v.dataFilter = dft;
v.dataFilter.subjectGroup  = group_pc1;
v.dataFilter.stims_allowed = 'all stimuli';
v.dataFilter.trial_sorting = 'chronological';
v.dataFilter.interval      = [.5, 20];
v.dataFilter.repetitions   = 1:5;

[~, oe_pc1, labs_pc1] = ModeSelector(v).extract;
L_pc1 = height(oe_pc1{1});
v.dataFilter.interval = [-22, -2];
[~, be_pc1] = ModeSelector(v).extract;
be_pc1 = cellfun(@(x) repmat(mean(x,1,'omitmissing'), L_pc1, 1, 1), be_pc1, 'UniformOutput', false);
% ev_pc1 = cellfun(@(b,o) o - b, be_pc1, oe_pc1, 'UniformOutput', false);
ev_pc1 = oe_pc1;

proj_pc1 = computeLDE(ev_pc1, labs_pc1, 'pooldata', true, 'nans2zeros', true, 'method', 'pca');

% mean PC1 per trial from pooled embedding (for scatter)
mean_pc1_all = squeeze(mean(proj_pc1.embedding{1}.reduction(:,1,:), 1, 'omitmissing'));
mean_pc1_all = mean_pc1_all(:);   % [nTrials_total x 1]

% per-subject embeddings for mean ± SEM, sign-aligned to pooled PC1
n_per_subj = cellfun(@(x) size(x,3), ev_pc1);
nSubj_pc1  = numel(ev_pc1);
pc1_mat = nan(max(n_per_subj),nSubj_pc1);
pc1_cols   = cell(1, nSubj_pc1);
ptr = 0;
for i = 1:nSubj_pc1
    proj_i     = computeLDE(ev_pc1(i), labs_pc1(i), 'pooldata', false, 'nans2zeros', true, 'method', 'pca');
    pc1_i      = squeeze(mean(proj_i.embedding{1}.reduction(:,1,:), 1, 'omitmissing'));
    pc1_i      = pc1_i(:);   % [ni x 1]
    pc1_mat(1:numel(pc1_i),i) = pc1_i;
end
nT_pc1  = size(pc1_mat, 1);

t_ax_pc1 = (1:nT_pc1)';
mu_pc1   = mean(pc1_mat, 2, 'omitmissing');
sem_pc1  = std(pc1_mat, [], 2, 'omitmissing') ./ sqrt(sum(isfinite(pc1_mat), 2));

hf_pc1 = figure; hold on;
scatter(repmat(t_ax_pc1, 1, nSubj_pc1), mean_pc1_all, 15, [.6 .6 .6], 'filled', 'MarkerFaceAlpha', 0.3);
errorbar(t_ax_pc1, mu_pc1, sem_pc1, '-', 'LineWidth', 1.5, 'CapSize', 0, 'Color', cfg.axcol);
xlabel('Trial number');  ylabel('Mean PC1');
title(group_pc1);  box off;
set(gca, 'color', cfg.bgcol, 'XColor', cfg.axcol, 'YColor', cfg.axcol);
set(gcf, 'color', cfg.bgcol);
cfg.figSize = 'small';  cfg.aspRatioType = 'square';  cfg.setLines = false;
cfg.setFigure;
cfg.saveFigure(hf_pc1, ['PC1 vs trial num ', group_pc1, ' raw'], saveType)

% R² between each subject's PC1 trajectory and trial number
t_vec = (1:nT_pc1)';
r2_subj = nan(nSubj_pc1, 1);
for i = 1:nSubj_pc1
    y  = pc1_mat(:, i);
    ok = isfinite(y);
    if sum(ok) >= 2
        r2_subj(i) = corr(t_vec(ok), y(ok))^2;
    end
end

figure;
histogram(r2_subj, 'FaceColor', cfg.axcol, 'EdgeAlpha', 0);
xlabel('R² (PC1 vs trial number)');  ylabel('Subjects');
box off;
set(gca, 'color', cfg.bgcol, 'XColor', cfg.axcol, 'YColor', cfg.axcol);
set(gcf, 'color', cfg.bgcol);
cfg.figSize = 'small';  cfg.aspRatioType = 'square';  cfg.setLines = false;
cfg.setFigure;
cfg.saveFigure(gcf, ['PC1 vs trial R2 hist ', group_pc1, ' raw'], saveType)

mean(r2_subj)

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