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
 

%% FIGURE 6

% stimulus dPCs except #1
v.dataFilter = dft;
v.dataFilter.mode_name = 'dpca';
v.dataFilter.mode_OI = 'stimulus';
v.dataFilter.mode_method = 'subtract';
v.dataFilter.subjectGroup = 'naïve';
v.dataFilter.mode_file = 'dpca_naive.mat';
hf = figure;
outMat_dn = v.plotDistancesHead;
v.dataFilter.mode_method = 'isolate';
hf = figure;
outMat_up = v.plotDistancesHead;
%
v.dataFilter = dft;

% stimulus dPCs except #1 (lower triangle = subtract, upper triangle = isolate)
v.dataFilter = dft;
v.dataFilter.mode_name = 'dpca';
v.dataFilter.mode_OI = 'stimulus';
v.dataFilter.mode_method = 'subtract';
v.dataFilter.subjectGroup = 'naïve';
v.dataFilter.mode_file = 'dpca_naive.mat';
hf = figure;
outMat_dn = v.plotDistancesHead;
close(hf)
v.dataFilter.mode_method = 'isolate';
hf = figure;
outMat_up = v.plotDistancesHead;
close(hf)
inMat = mean(outMat_dn.distMat3d,3,"omitmissing");
outMat_up = mean(outMat_up.distMat3d,3,"omitmissing");
triu_idx = logical(triu(true(size(outMat_up))));
inMat(triu_idx) = outMat_up(triu_idx);
%plot
labs = outMat_dn.all_labs{1};
ntrials = numel(labs);
hf = figure;
imagesc(1-inMat)
axis square; hold on
xticks(1:ntrials); xticklabels(labs); xtickangle(90)
yticks(1:ntrials); yticklabels(labs)
xlabel('Stimulus type')
ylabel('Stimulus type')
clim([.2 1])
colormap(cfg.colormapName)
a = colorbar('Color',cfg.axcol);
a.Label.String = 'correlation';
a.Label.FontSize= gca().FontSize;
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
set(gcf, 'color', cfg.bgcol); 
hold off
%
v.dataFilter = dft;


% only novelty dPC
%
% reconstructed data iFR over experiment time (imagesc)
v.dataFilter = dft;
v.dataFilter.mode_name = 'dpca';
v.dataFilter.mode_OI = 'novelty';
v.dataFilter.mode_method = 'isolate';
v.dataFilter.interval = [];
v.dataFilter.subjectGroup = 'naïve';
v.dataFilter.mode_file = 'dpca_naive.mat';
closeup_interval = [-1 3];
[~,events,labs] = ModeSelector(v).extract;
events = cellfun(@(x) squeeze(mean(x,2,'omitmissing')),events,'UniformOutput',false);
M = [];
for i = 1:numel(events)
M(:,:,i) = events{i};
end
inMat = mean(M,3,'omitmissing');
%plot
labs = outMat_dn.all_labs{1};
ntrials = numel(labs);
hf1 = figure;
y = 1:ntrials;
fs = v.filtered_traces{1}.framerate;
x = [0:1/fs:height(inMat)/fs-1/fs] -30;
imagesc(x,y,inMat')
yticks(1:ntrials); yticklabels(labs);
ylabel('Stimulus type')
xlabel('Time [s]')
% clim([0 1])
colormap(cfg.colormapName)
a = colorbar('Color',cfg.axcol);
a.Label.String = 'iFR [Hz]';
a.Label.FontSize= gca().FontSize;
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
set(gcf, 'color', cfg.bgcol); 
% timecourse
hf2 = figure;
y = mean(inMat,2,"omitmissing");
plot(x,y,'Color',cfg.c(1,:),'LineWidth',1)
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
set(gcf, 'color', cfg.bgcol); 
ylabel('iFR [Hz]')
xlabel('Time [s]')
axis tight, box off
ylim([0 max(ylim)])
% heatmap close-up
figure(hf1);
xlim(closeup_interval)
% timecourse over reps
v.dataFilter.interval = t_lim_sec;
hf2 = figure;
for i = 1:5
    v.dataFilter.repetitions = i;
    [~,events,labs] = ModeSelector(v).extract;
    events = cellfun(@(x) squeeze(mean(x,2,'omitmissing')),events,'UniformOutput',false);
    M = [];
    for j = 1:numel(events)
    M(:,:,j) = events{j};
    end
    inMat = mean(M,3,'omitmissing');
    y = mean(inMat,2,"omitmissing");
    x = [0:1/fs:height(inMat)/fs-1/fs] +min(v.dataFilter.interval);
    plot(x,y,'Color',cfg.c(i,:),'LineWidth',2); hold on
    set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
    set(gcf, 'color', cfg.bgcol); 
    ylabel('iFR [Hz]')
    xlabel('Time [s]')
    axis tight, box off
    ylim([0 max(ylim)])
end
legend({'rep1','rep2','rep3','rep4','rep5'})
%
v.dataFilter = dft;


% unit firing distributions (recontructed data from novelty mode)
%
yrange = [0 .05];
v.dataFilter = dft;
v.dataFilter.mode_name = 'dpca';
v.dataFilter.mode_OI = 'novelty';
v.dataFilter.mode_method = 'isolate';
v.dataFilter.interval = [0 3];
% naive
v.dataFilter.subjectGroup = 'naïve';
v.dataFilter.mode_file = 'dpca_naive.mat';
v.dataFilter.stims_allowed = 'all stimuli';
hf = figure;
c = v.plotTrialActivityMetricHead('method','avg intensity', 'plotType', 'boxplot_repetitions');
ylim(yrange)
[~,~,stats] = kruskalwallis(c);
c = multcompare(stats)
% "trained stimuli"
v.dataFilter.stims_allowed = {'Arg','Ala','His'};
hf = figure;
v.plotTrialActivityMetricHead('method','avg intensity', 'plotType', 'boxplot_repetitions');
ylim(yrange)
% "novel stimuli"
v.dataFilter.stims_allowed = {'Trp','Ser','Leu'};
hf = figure;
v.plotTrialActivityMetricHead('method','avg intensity', 'plotType', 'boxplot_repetitions');
ylim(yrange)
%
% trained
v.dataFilter.subjectGroup = 'trained';
v.dataFilter.mode_file = 'dpca_trained.mat';
v.dataFilter.stims_allowed = 'all stimuli';
hf = figure;
v.plotTrialActivityMetricHead('method','avg intensity', 'plotType', 'boxplot_repetitions');
ylim(yrange)
% "trained stimuli"
v.dataFilter.stims_allowed = {'Arg','Ala','His'};
hf = figure;
v.plotTrialActivityMetricHead('method','avg intensity', 'plotType', 'boxplot_repetitions');
ylim(yrange)
% "novel stimuli"
v.dataFilter.stims_allowed = {'Trp','Ser','Leu'};
hf = figure;
v.plotTrialActivityMetricHead('method','avg intensity', 'plotType', 'boxplot_repetitions');
ylim(yrange)
%
% uncoupled
v.dataFilter.subjectGroup = 'uncoupled';
v.dataFilter.mode_file = 'dpca_uncoupled.mat';
v.dataFilter.stims_allowed = 'all stimuli';
hf = figure;
v.plotTrialActivityMetricHead('method','avg intensity', 'plotType', 'boxplot_repetitions');
ylim(yrange)
% "trained stimuli"
v.dataFilter.stims_allowed = {'Arg','Ala','His'};
hf = figure;
v.plotTrialActivityMetricHead('method','avg intensity', 'plotType', 'boxplot_repetitions');
ylim(yrange)
% "novel stimuli"
v.dataFilter.stims_allowed = {'Trp','Ser','Leu'};
hf = figure;
v.plotTrialActivityMetricHead('method','avg intensity', 'plotType', 'boxplot_repetitions');
ylim(yrange)


% intertrial correlations (units vs stimulus dPCs)
%
v.dataFilter = dft;
metric = 'correlation';
v.dataFilter.subjectGroup = 'naïve';
[~, units3d, units_vals] = plotTuningCorrelationsOverReps(v, metric); % native units
v.dataFilter.mode_name = 'dpca';
v.dataFilter.mode_OI = 'stimulus';
v.dataFilter.mode_method = 'remove';
v.dataFilter.mode_file = 'dpca_naive.mat';
[~, dPC3d, dPC_vals] = plotTuningCorrelationsOverReps(v, metric); % stimulus dPCs
% multicolor histogram
hf = figure;
edges = -1:.05:1;
histogram(units_vals,edges,'FaceColor',cfg.c(1,:))
hold on
histogram(dPC_vals,edges,'FaceColor',cfg.c(2,:))
line([.5,.5],[0 300],'Color','r','LineStyle','--','LineWidth',2)
axis square
box off
xlabel(['avg inter-rep ',metric])
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
set(gcf, 'color', cfg.bgcol); 
[~,p] = kstest2(units_vals,dPC_vals);
disp(['2-sample KS test for unit vs dPC intertrial corrs - pval: ',num2str(p)])


