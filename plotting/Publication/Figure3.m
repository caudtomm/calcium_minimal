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
cfg = PlotConfig('colormapName','lapaz','favouriteColors',[84,85,73]); % (test1, test2, ctrl)
cfg.savePath = savepath;

v = ExperimentViewer(experiment);
v.plotConfig = cfg;
cfg = v.plotConfig;

v.dataFilter.traceType = 'pSpike';
v.dataFilter.interval = [1,20];
v.dataFilter.trial_sorting = 'stim_id';
dft = v.dataFilter;
 

%% FIGURE 3

% general suppression score
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'naïve';
c = v.plotUnitActivityMetricHead('method','general suppression score');
% histogram
figure; histogram(cell2mat(c),100, 'FaceColor','k','EdgeAlpha',0);
box off; axis square
xlabel('suppression [a.u]'); ylabel('histogram')
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
set(gcf, 'color', cfg.bgcol);

% selectivity of tuning
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

% drift strength vs activity
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'naïve';
R = driftMetricsFigure(v,'plot',false);
drift = R.driftStrength3d;
% vs stimulus-specific activity on the first trial of the pair
hf = figure; 
v.dataFilter.repetitions = 1;
out = v.plotUnitActivityMetricHead('metric','avg intensity');
firing = cell2mat(out);
thisdrift = drift(:,:,1);
subplot(211); plotHeatmapAndIsoclines(firing(:),thisdrift(:),50,1,0,1);
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
subplot(212); plotHeatmapAndIsoclines(firing(:),thisdrift(:),50,1,0,1);
xlim([0 1])
ylim([-15 15])
colormap(cfg.colormapName)
title('rep 4->5')
xlabel('iFR on rep 4 [Hz]'); ylabel('drift modulus')
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
% vs averaged activity across stimuli and reps
hf = figure;
v.dataFilter.repetitions = 1:5;
out = v.plotUnitActivityMetricHead('metric','avg intensity');
firing = cell2mat(out); firing = mean(firing,[2,3],'omitmissing');
thisdrift = drift(:,:,1); thisdrift = mean(thisdrift,2,'omitmissing');
subplot(211); plotHeatmapAndIsoclines(firing(:),thisdrift(:),50,1,0,1);
xlim([0 .4])
ylim([-2.5 2.5])
colormap(cfg.colormapName)
title('rep 1->2')
xlabel('average iFR [Hz]'); ylabel('average drift modulus')
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
thisdrift = drift(:,:,4); thisdrift = mean(thisdrift,2,'omitmissing');
subplot(212); plotHeatmapAndIsoclines(firing(:),thisdrift(:),50,1,0,1);
xlim([0 .4])
ylim([-2.5 2.5])
colormap(cfg.colormapName)
title('rep 4->5')
xlabel('average iFR [Hz]'); ylabel('average drift modulus')
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);

% stability of tuning
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


% unit firing distributions
yrange = [0 .2];
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'naïve';
v.dataFilter.trial_sorting = 'chronological';
% over trials
v.dataFilter.stims_allowed = 'all trials';
hf = figure;
v.plotTrialActivityMetricHead('method','avg intensity')
ylim(yrange)
% over repetitions
v.dataFilter.stims_allowed = 'all stimuli';
hf = figure;
v.plotTrialActivityMetricHead('method','avg intensity', 'plotType', 'boxplot_repetitions')
ylim(yrange)
% over repetitions (Leu)
v.dataFilter.stims_allowed = {'Leu'};
hf = figure;
v.plotTrialActivityMetricHead('method','avg intensity', 'plotType', 'boxplot_repetitions')
ylim(yrange)

% unit firing distributions (baseline)
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'naïve';
v.dataFilter.trial_sorting = 'chronological';
v.dataFilter.interval = [-22 -2];
% over trials
v.dataFilter.stims_allowed = 'all trials';
hf = figure;
v.plotTrialActivityMetricHead('method','avg intensity')
ylim(yrange)
% over repetitions
v.dataFilter.stims_allowed = 'all stimuli';
hf = figure;
v.plotTrialActivityMetricHead('method','avg intensity', 'plotType', 'boxplot_repetitions')
ylim(yrange)
% over repetitions (Leu)
v.dataFilter.stims_allowed = {'Leu'};
hf = figure;
v.plotTrialActivityMetricHead('method','avg intensity', 'plotType', 'boxplot_repetitions')
ylim(yrange)

