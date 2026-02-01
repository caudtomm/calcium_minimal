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

% odor manifolds (all odors)
folder_tag = 'odors'; hf = plotMetrics(v, folder_tag, cfg);
savePlot(hf(1),'Capacity',[0 .15], folder_tag, cfg, saveType);
savePlot(hf(2),'Dimension',[0 30], folder_tag, cfg, saveType);
savePlot(hf(3),'Radius',[0 2], folder_tag, cfg, saveType);
savePlot(hf(5),'Center alignment', folder_tag,[0 1], cfg, saveType);
savePlot(hf(6),'Axis alignment', folder_tag,[0 .5], cfg, saveType);

% repetition manifolds
folder_tag = 'repetitions'; hf = plotMetrics(v, folder_tag, cfg);
savePlot(hf(1),'Capacity',[0 .15], folder_tag, cfg, saveType);





function savePlot(hf, metric_name, yrange, folder_tag, cfg, saveType)
    figure(hf);
    ylim(yrange)
    title(''); ylabel(metric_name)
    cfg.figSize = 'tiny';
    cfg.aspRatioType = 'tall';
    cfg.lineWidth = .5;
    cfg.setLines = true;
    cfg.setFigure;
    cfg.saveFigure(gcf,['gcmc ',folder_tag, ' ',metric_name,' boxplot'], saveType)
end

function hf = plotMetrics(v, folder_tag, cfg)
    indir = fullfiletol('manifold_data',folder_tag);
    all_results = GCMC_Analysis(v).extractResultsFromMultipleSubjects(fullfiletol(indir,'odorexp004_IC1_130625'));
    group_data = GCMC_Analysis(v).clusterByGroup(all_results);
    new_group_data = mergeTrainedGroups(group_data);
    
    hf = GCMC_Plotting.plotBoxplotsByGroup(new_group_data, cfg, false); 
end

function new_group_data = mergeTrainedGroups(group_data)
    % merge trained groups
    new_group_data(1) = group_data(1);
    new_group_data(2) = group_data(5);
    new_group_data(3).group_name = 'trained';
    new_group_data(3).data = group_data(2).data;
    new_group_data(3).data = [new_group_data(3).data; group_data(3).data];
    new_group_data(3).data = [new_group_data(3).data; group_data(4).data];
end


%% template matching lines

v.dataFilter = dft;
grouptag = 'uncoupled';
stimtag = 'allstims';
v.dataFilter.subjectGroup = 'uncoupled';
v.dataFilter.stims_allowed = 'all stimuli';
v.dataFilter.interval = [.5, 20];
v.dataFilter.repetitions = 2:5;
method = 'correlation';
focus_stims = 'all stimuli';

hf = figure;
c = v.plotDiscriminationHead(...
    'plotType', 'performance_lines', ...
    'method',method, ...
    'focus_stims', focus_stims, ...
    'zscore',false);
for i = 2:width(c)
    p = signrank(c(:,1),c(:,i));
    disp(['Paired Wilcoxon signed-rank test - reps 1 vs ', num2str(i),': ',num2str(p)])
end
xlabel('Template trial #')
ylabel('Performance')
ylim([0 1])
cfg.figSize = "tiny";
cfg.aspRatioType = "tall";
cfg.setFigure;
cfg.saveFigure(gcf,['template ',grouptag, ' ',stimtag,' lines 2-5'], saveType)


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