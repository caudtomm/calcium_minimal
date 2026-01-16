dbstop if error

s = true; % save figures to files?
savepath = 'bin3';
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
cfg.savePath = savepath;

v = ExperimentViewer(experiment);
v.plotConfig = cfg;
cfg = v.plotConfig;

v.dataFilter.traceType = 'pSpike';
v.dataFilter.interval = [1,20];
v.dataFilter.trial_sorting = 'stim_id';
dft = v.dataFilter;
 

%% FIGURE 3


%% unit firing distributions
baseline_interval = [-22 -2];
odor_interval = dft.interval;
yrange = [-.2 .3];
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'trained';
v.dataFilter.trial_sorting = 'chronological';
% over trials
v.dataFilter.stims_allowed = 'all trials';
out = v.plotUnitActivityMetricHead('method','avg intensity');
v.dataFilter.interval = baseline_interval;
outbase = v.plotUnitActivityMetricHead('method','avg intensity');
data = cell2mat(out) - cell2mat(outbase);
[~,~, labs] = ModeSelector(v).extract; labs = labs{1};
hf = figure; b = prettyBoxplot(data,labs,'scatterSize',5,'plotLine',true);
ylabel('Cellwise delta iFR (Hz)')
ylim(yrange)
cfg.figSize = "large";
cfg.aspRatioType = "wide";
cfg.lineWidth = .5;
cfg.setFigure
cfg.saveFigure(gcf,'trained unit delta firing distribution', saveType)
cfg = v.plotConfig;

% unit firing distributions (repetitions)
baseline_interval = [-22 -2];
odor_interval = dft.interval;
yrange = [-.2 .3];
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'trained';
v.dataFilter.trial_sorting = 'chronological';
v.dataFilter.stims_allowed = 'all stimuli';
[~,~, labs] = ModeSelector(v).extract; labs = labs{1};
% over reps
stims = unique(labs); nstims = numel(stims);
data = [];
for i = 1:nstims
    v.dataFilter.stims_allowed = stims(i);
    v.dataFilter.interval = odor_interval;
    out = v.plotUnitActivityMetricHead('method','avg intensity');
    v.dataFilter.interval = baseline_interval;
    outbase = v.plotUnitActivityMetricHead('method','avg intensity');
    data = [data; cell2mat(out) - cell2mat(outbase)];
end
hf = figure;
subplot(121); b = prettyBoxplot(data,{'1','2','3','4','5'},'scatterSize',5,'plotLine',true);
ylabel('Cellwise delta iFR (Hz)')
ylim(yrange)
% Leu
v.dataFilter.stims_allowed = {'Leu'};
v.dataFilter.interval = odor_interval;
out = v.plotUnitActivityMetricHead('method','avg intensity');
v.dataFilter.interval = baseline_interval;
outbase = v.plotUnitActivityMetricHead('method','avg intensity');
data = cell2mat(out) - cell2mat(outbase);
subplot(122); b = prettyBoxplot(data,{'1','2','3','4','5'},'scatterSize',5,'plotLine',true);
ylim(yrange)
ylabel('Cellwise delta iFR (Hz)')
cfg.figSize = "medium";
cfg.aspRatioType = "square";
cfg.lineWidth = .5;
cfg.setFigure
cfg.saveFigure(gcf,'trained unit delta firing distribution repetitions', saveType)
cfg = v.plotConfig;

% 
% % unit firing distributions
% baseline_interval = [-22 -2];
% odor_interval = dft.interval;
% yrange = [0 .2];
% v.dataFilter = dft;
% v.dataFilter.subjectGroup = 'naïve';
% v.dataFilter.trial_sorting = 'chronological';
% % over trials
% v.dataFilter.stims_allowed = 'all trials';
% hf = figure;
% out = v.plotTrialActivityMetricHead('method','avg intensity');
% ylim(yrange)
% % over repetitions
% v.dataFilter.stims_allowed = 'all stimuli';
% hf = figure;
% v.plotTrialActivityMetricHead('method','avg intensity', 'plotType', 'boxplot_repetitions')
% ylim(yrange)
% % over repetitions (Leu)
% v.dataFilter.stims_allowed = {'Leu'};
% hf = figure;
% v.plotTrialActivityMetricHead('method','avg intensity', 'plotType', 'boxplot_repetitions')
% ylim(yrange)
% 

%% general suppression score histogram
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'naïve';
c = v.plotUnitActivityMetricHead('method','general suppression score');
% histogram
figure; histogram(cell2mat(c),100, 'FaceColor','k','EdgeAlpha',0);
box off; axis square
xlabel('suppression [a.u]'); ylabel('histogram')
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
set(gcf, 'color', cfg.bgcol);

%% selectivity of tuning histogram
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'naïve';
c = v.plotUnitActivityMetricHead('method','selectivity of tuning');
% histogram
figure; histogram(cell2mat(c),100, 'FaceColor','k','EdgeAlpha',0);
box off; axis square
xlim([0 1])
xlabel('tuning selectivity [a.u]'); ylabel('histogram')
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
set(gcf, 'color', cfg.bgcol);

%% stability of tuning histogram
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'naïve';
c = v.plotUnitActivityMetricHead('method','stability of tuning');
% histogram
figure; histogram(cell2mat(c),100, 'FaceColor','k','EdgeAlpha',0);
box off; axis square
xlim([-1 1])
xlabel('tuning stability [a.u]'); ylabel('histogram')
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
set(gcf, 'color', cfg.bgcol);


%% drift strength vs activity
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'naïve';
R = driftMetricsFigure(v,'plot',true);
drift = R.driftStrength3d;
driftAx = R.figures.drift_cdf.Axes(1);
hf = figure;
% vs averaged activity across stimuli and reps
v.dataFilter.repetitions = 1:5;
out = v.plotUnitActivityMetricHead('metric','avg intensity');
firing = cell2mat(out); firing = mean(firing,[2,3],'omitmissing');
thisdrift = drift(:,:,1); thisdrift = mean(thisdrift,2,'omitmissing');
subplot(323); plotHeatmapAndIsoclines(firing(:),thisdrift(:),50,1,0,1);
xlim([0 .4])
ylim([-2.5 2.5])
colormap(cfg.colormapName)
title('rep 1->2')
xlabel('average iFR [Hz]'); ylabel('average drift modulus')
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
thisdrift = drift(:,:,4); thisdrift = mean(thisdrift,2,'omitmissing');
subplot(324); plotHeatmapAndIsoclines(firing(:),thisdrift(:),50,1,0,1);
xlim([0 .4])
ylim([-2.5 2.5])
colormap(cfg.colormapName)
title('rep 4->5')
xlabel('average iFR [Hz]'); ylabel('average drift modulus')
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
% vs stimulus-specific activity on the first trial of the pair
v.dataFilter.repetitions = 1;
out = v.plotUnitActivityMetricHead('metric','avg intensity');
firing = cell2mat(out);
thisdrift = drift(:,:,1);
subplot(325); plotHeatmapAndIsoclines(firing(:),thisdrift(:),50,1,0,1);
xlim([0 1])
ylim([-15 15])
colormap(cfg.colormapName)
title('rep 1->2')
xlabel('iFR on rep 1 [Hz]'); ylabel('drift modulus')
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
v.dataFilter.repetitions = 4;
out = v.plotUnitActivityMetricHead('metric','avg intensity');
firing = cell2mat(out);
thisdrift = drift(:,:,v.dataFilter.repetitions);
subplot(326); plotHeatmapAndIsoclines(firing(:),thisdrift(:),50,1,0,1);
xlim([0 1])
ylim([-15 15])
colormap(cfg.colormapName)
title('rep 4->5')
xlabel('iFR on rep 4 [Hz]'); ylabel('drift modulus')
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);

cfg.figSize = "large";
cfg.aspRatioType = "square";
cfg.setFigure;
cfg.saveFigure(gcf,'naive drift vs activity', saveType)
cfg = v.plotConfig;

