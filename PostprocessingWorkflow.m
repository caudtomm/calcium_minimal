%% Knobs

s = false; % save figures to files?

% for sliding windows
window_duration = 1; % [seconds]
t_lim_sec = [-5 35]; % from 5 sec before to 35 seconds after stimulus onset
overlap = .3; % [seconds]

%% Load dataset
filename = 'odorexp004_IC1_130625.mat';

%
experiment = load(filename).a; % Experiment object
% or
experiment = a; clear a

% avoid any spelling mismatches
for i = 1:numel(experiment.traces)
experiment.traces{i}.subject_group = experiment.subjectTab.group{i};
end


%% initialize output figure saving
cfg = PlotConfig('theme', 'light');

figs = FigureSaver;
figs.outputfolder = fullfiletol('figures',extractBefore(filename,'.'));
figs.outputfile = 'all_plots.pdf';
figs.config = cfg;

%% Plotting average similarity matrices and related metrics for each experimental group.
v = ExperimentViewer(experiment);
v.dataFilter.traceType = 'pSpike';
v.plotConfig = cfg;

% plot for naive fish
v.dataFilter.subjectGroup = 'naïve';
hf = figure;
v.plotDistances('ps_lim',[1,20],'trial_sorting','stim_id');
figs.title = 'Naive Group';
figs.append(hf);
close(hf)

% plot for trained fish
v.dataFilter.subjectGroup = 'trained';
hf = figure;
v.plotDistances('ps_lim',[1,20],'trial_sorting','stim_id');
figs.title = 'Trained Groups';
figs.append(hf);
close(hf)

% plot for uncoupled fish
v.dataFilter.subjectGroup = 'uncoupled';
hf = figure;
v.plotDistances('ps_lim',[1,20],'trial_sorting','stim_id');
figs.title = 'Uncoupled Group';
figs.append(hf);
close(hf)

% stimulus repetition comparisons
hf = plotRepetitionDistances(v,[1 20],'correlation'); % outputs 2 figures
figs.title = 'Repetition comparison';
figs.append(hf);
close(hf)

% stimulus repetitions, second by second
windows = defineTimeWindows(window_duration,t_lim_sec,overlap);
[hf,data] = similarityDynamics(v,figs,windows,s);
figs.title = 'Repetitions over time';
figs.append(hf);
close(hf)

%% Discrimination analysis (template-matching)

% template matching, stimulus window
hf = plotDiscriminationPerformanceMats(v,[1 20], 'correlation','all trials',1:5,false); % outputs 1 figure
figs.title = 'Template-match performance comparison';
figs.append(hf);
close(hf)

% template matching, second by second
windows = defineTimeWindows(window_duration,t_lim_sec,overlap);
[hf,data] = discriminationDynamics(v,figs,windows,s);
figs.title = 'Discrimination over time';
figs.append(hf);
close(hf)

% template matching, stimulus window, focus on performance for novel
% stimuli
hf = plotDiscriminationPerformanceMats(v,[1 20], 'correlation','all novel',1:5,false); % outputs 1 figure
figs.title = 'Template-match performance on novel stimuli';
figs.append(hf);
close(hf)

% template matching, stimulus window, focus on performance for Leu
hf = plotDiscriminationPerformanceMats(v,[1 20], 'correlation',{'Leu'},1:5,false); % outputs 1 figure
figs.title = 'Template-match performance on Leu';
figs.append(hf);
close(hf)

% template matching, stimulus window, focus on performance for familiar
% stimuli
hf = plotDiscriminationPerformanceMats(v,[1 20], 'correlation','all familiar',1:5,false); % outputs 1 figure
figs.title = 'Template-match performance on familiar stimuli';
figs.append(hf);
close(hf)

% template matching, stimulus window, focus on performance for CS+
% stimuli
hf = plotDiscriminationPerformanceMats(v,[1 20], 'correlation','all CS+',1:5,false); % outputs 1 figure
figs.title = 'Template-match performance on CS+ stimuli';
figs.append(hf);
close(hf)

% template matching, stimulus window, focus on performance for CS-
% stimuli
hf = plotDiscriminationPerformanceMats(v,[1 20], 'correlation','all CS-',1:5,false); % outputs 1 figure
figs.title = 'Template-match performance on CS- stimuli';
figs.append(hf);
close(hf)

% template matching based only on trials 1:4, stimulus window
hf = plotDiscriminationPerformanceMats(v,[1 20], 'correlation','all trials',1:4,false); % outputs 1 figure
figs.title = 'Template-match performance comparison';
figs.append(hf);
close(hf)

% template matching based only on trials 2:4, stimulus window
hf = plotDiscriminationPerformanceMats(v,[1 20], 'correlation','all trials',2:4,false); % outputs 1 figure
figs.title = 'Template-match performance comparison';
figs.append(hf);
close(hf)

% template matching based only on trials 2:5, stimulus window
hf = plotDiscriminationPerformanceMats(v,[1 20], 'correlation','all trials',2:5,false); % outputs 1 figure
figs.title = 'Template-match performance comparison';
figs.append(hf);
close(hf)

%% 

[hf, out, groups, odor_sets] = plotTrialMetricFigure(v, [1 20], 'normalized population sparseness');

[hf, out, groups, odor_sets] = plotTrialMetricFigure(v, [1 20], 'participation ratio');

[hf, out, groups, odor_sets] = plotTrialMetricFigure(v, [1 20], 'max intensity', 'cells');

[hf, out, groups, odor_sets] = plotTrialMetricFigure(v, [1 20], 'avg intensity', 'cells');

[hf, out, groups, odor_sets] = plotTrialMetricFigure(v, [1 20], 'variance', 'cells');

[hf, out, groups, odor_sets] = plotTrialMetricFigure(v, [1 20], 'max intensity', 'frames');

[hf, out, groups, odor_sets] = plotTrialMetricFigure(v, [1 20], 'avg intensity', 'frames');

[hf, out, groups, odor_sets] = plotTrialMetricFigure(v, [1 20], 'variance', 'frames');

%%

plotUnitTuningFigure(v, [1 20], 'pca')
plotUnitTuningFigure(v, [1 3], 'pca')
plotUnitTuningFigure(v, [3 8], 'pca')
plotUnitTuningFigure(v, [12 20], 'pca')

plotUnitTuningFigure(v, [1 20], 'isomap')
plotUnitTuningFigure(v, [1 3], 'isomap')
plotUnitTuningFigure(v, [3 8], 'isomap')
plotUnitTuningFigure(v, [12 20], 'isomap')

%% 

% Stability of Tuning depends on correlation- or cosine-similarity-based
% definition (see extractActivityMetric.m)

v.dataFilter.subjectGroup = 'naïve';
compareModeMetricsFigure(v);
v.dataFilter.subjectGroup = 'trained';
compareModeMetricsFigure(v);
v.dataFilter.subjectGroup = 'uncoupled';
compareModeMetricsFigure(v);

%%






