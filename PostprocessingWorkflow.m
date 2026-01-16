%% Knobs

dbstop if error

s = false; % save figures to files?

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
cfg = PlotConfig('theme', 'light');

figs = FigureSaver;
figs.outputfolder = fullfiletol('figures',extractBefore(filename,'.'));
figs.outputfile = 'figure2.pdf';
figs.config = cfg;


v = ExperimentViewer(experiment);
v.plotConfig = cfg;

v.dataFilter.traceType = 'pSpike';
v.dataFilter.interval = [1,20];
v.dataFilter.trial_sorting = 'stim_id';
dft = v.dataFilter;

%% GCMC: extract manifolds and save to .mat files for Python analysis

outdir = fullfiletol('manifold_data','trial_slide_windows', extractBefore(filename,'.'));

% filter data
v.dataFilter.traceType = 'pSpike';
v.dataFilter.subjectIDs = {};
v.dataFilter.subjectGroup = 'all';
v.dataFilter.interval = [.5 25];
v.dataFilter.repetitions = [1:5];
v.dataFilter.stims_allowed = 'all stimuli';
v.dataFilter.mode_name = 'native_units';
v.dataFilter.mode_method = 'mode_values';
v.dataFilter.mode_OI = 'all';
v.dataFilter.mode_file = '';

% GCMC_Analysis(v).outputDataFiles('odor',outdir); % each manifold: one stimulus, all repetitions, one subject

[pathlist, manifolds] = GCMC_Analysis(v).outputDataFiles_SlidingWindow('trial',outdir,[-6,36],4,2); % each manifold: one stimulus, all repetitions, one subject

v.dataFilter = dft; % recover

%% GCMC: load capacity results from .mat files for MATLAB analysis

indir = fullfiletol('manifold_data', extractBefore(filename,'.'), 'results');
indir = fullfiletol('manifold_data','odors_man_windows'); % tempoarily use this folder for the test dataset
[~, avg_results] = GCMC_Analysis(v).extractResults(indir); % # TODO this doesn't take into account multiple subjects yet
% # TODO save results to ActivityTraces inside v
% # TODO some plotting here

GCMC_Analysis.plotBoxplotsForEachMetric(avg_results,cfg);
GCMC_Analysis.plotMetricStability(results, cfg)

all_results = GCMC_Analysis(v).extractResultsFromMultipleSubjects(indir);
group_data = GCMC_Analysis(v).clusterByGroup(all_results);

% plotting
GCMC_Analysis.plotBoxplotsByGroup(new_group_data, cfg, false)

%
all_results = GCMC_Analysis(v).extractResults_SlidingWindow(indir,'odorexp004_IC1_130625');

%% Plotting average similarity matrices and related metrics for each experimental group.

% plot for naive fish
v.dataFilter.subjectGroup = 'naïve';
hf = figure;
v.plotDistancesHead;
v.dataFilter = dft;
figs.title = 'Naive Group';
figs.append(hf);
close(hf)

% plot for trained fish
v.dataFilter.subjectGroup = 'trained';
hf = figure;
v.plotDistancesHead;
v.dataFilter = dft;
figs.title = 'Trained Groups';
figs.append(hf);
close(hf)

% plot for uncoupled fish
v.dataFilter.subjectGroup = 'uncoupled';
hf = figure;
v.plotDistancesHead;
v.dataFilter = dft;
figs.title = 'Uncoupled Group';
figs.append(hf);
close(hf)

% stimulus repetition comparisons
hf = plotRepetitionDistances(v,'correlation'); % outputs 2 figures
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
hf = plotDiscriminationPerformanceMats(v, 'correlation','all trials',false); % outputs 1 figure
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
hf = plotDiscriminationPerformanceMats(v, 'correlation','all novel',false); % outputs 1 figure
figs.title = 'Template-match performance on novel stimuli';
figs.append(hf);
close(hf)

% template matching, stimulus window, focus on performance for Leu
hf = plotDiscriminationPerformanceMats(v, 'correlation',{'Leu'},false); % outputs 1 figure
figs.title = 'Template-match performance on Leu';
figs.append(hf);
close(hf)

% template matching, stimulus window, focus on performance for familiar
% stimuli
hf = plotDiscriminationPerformanceMats(v, 'correlation','all familiar',false); % outputs 1 figure
figs.title = 'Template-match performance on familiar stimuli';
figs.append(hf);
close(hf)

% template matching, stimulus window, focus on performance for CS+
% stimuli
hf = plotDiscriminationPerformanceMats(v, 'correlation','all CS+',false); % outputs 1 figure
figs.title = 'Template-match performance on CS+ stimuli';
figs.append(hf);
close(hf)

% template matching, stimulus window, focus on performance for CS-
% stimuli
hf = plotDiscriminationPerformanceMats(v, 'correlation','all CS-',false); % outputs 1 figure
figs.title = 'Template-match performance on CS- stimuli';
figs.append(hf);
close(hf)

% template matching based only on trials 1:4, stimulus window
v.dataFilter.repetitions = 1:4;
hf = plotDiscriminationPerformanceMats(v, 'correlation','all trials',false); % outputs 1 figure
v.dataFilter = dft;
figs.title = 'Template-match performance comparison';
figs.append(hf);
close(hf)

% template matching based only on trials 2:4, stimulus window
v.dataFilter.repetitions = 2:4;
hf = plotDiscriminationPerformanceMats(v, 'correlation','all trials',false); % outputs 1 figure
v.dataFilter = dft;
figs.title = 'Template-match performance comparison';
figs.append(hf);
close(hf)

% template matching based only on trials 2:5, stimulus window
v.dataFilter.repetitions = 2:5;
hf = plotDiscriminationPerformanceMats(v, 'correlation','all trials',false); % outputs 1 figure
v.dataFilter = dft;
figs.title = 'Template-match performance comparison';
figs.append(hf);
close(hf)

%% 

[hf, out, groups, odor_sets] = plotTrialMetricFigure(v, 'normalized population sparseness');

[hf, out, groups, odor_sets] = plotTrialMetricFigure(v, 'participation ratio');

[hf, out, groups, odor_sets] = plotTrialMetricFigure(v, 'max intensity', 'cells');

[hf, out, groups, odor_sets] = plotTrialMetricFigure(v, 'avg intensity', 'cells');

[hf, out, groups, odor_sets] = plotTrialMetricFigure(v, 'variance', 'cells');

[hf, out, groups, odor_sets] = plotTrialMetricFigure(v, 'max intensity', 'frames');

[hf, out, groups, odor_sets] = plotTrialMetricFigure(v, 'avg intensity', 'frames');

[hf, out, groups, odor_sets] = plotTrialMetricFigure(v, 'variance', 'frames');

%%

% embedding tuning curves using PCA
v.dataFilter.interval = [1 20];
[hf, out] = plotUnitTuningFigure(v,'pca');
v.dataFilter.interval = [1 3];
plotUnitTuningFigure(v, 'pca')
v.dataFilter.interval = [3 8];
plotUnitTuningFigure(v, 'pca')
v.dataFilter.interval = [12 20];
plotUnitTuningFigure(v, 'pca')

% embedding tuning curves using Isomap
v.dataFilter.interval = [1 20];
plotUnitTuningFigure(v, 'isomap')
v.dataFilter.interval = [1 3];
plotUnitTuningFigure(v, 'isomap')
v.dataFilter.interval = [3 8];
plotUnitTuningFigure(v, 'isomap')
v.dataFilter.interval = [12 20];
plotUnitTuningFigure(v, 'isomap')

v.dataFilter = dft; % recover

%% 

% Stability of Tuning depends on correlation- or cosine-similarity-based
% definition (see extractActivityMetric.m)

v.dataFilter.subjectGroup = 'naïve';
% v.dataFilter.mode_file = 'dpca_naive.mat';
compareModeMetricsFigure(v);

v.dataFilter.subjectGroup = 'trained';
% v.dataFilter.mode_file = 'dpca_trained.mat';
compareModeMetricsFigure(v);

v.dataFilter.subjectGroup = 'uncoupled';
% v.dataFilter.mode_file = 'dpca_uncoupled.mat';
compareModeMetricsFigure(v);

v.dataFilter = dft; % recover

%%

v.dataFilter = dft; % recover
v.dataFilter.mode_file = '';
v.dataFilter.mode_name = 'pca';
v.dataFilter.mode_OI = [1 2 3];
v.dataFilter.mode_method = 'mode_values';

v.dataFilter.subjectGroup = 'naïve';
driftMetricsFigure(v);

v.dataFilter.subjectGroup = 'trained';
driftMetricsFigure(v);

v.dataFilter.subjectGroup = 'uncoupled';
driftMetricsFigure(v);

v.dataFilter = dft; % recover

%% 

metric = 'correlation';


% native units
v.dataFilter.mode_file = '';
v.dataFilter.mode_name = 'native_units';
v.dataFilter.mode_OI = 'all';
v.dataFilter.mode_method = 'mode_values';

v.dataFilter.subjectGroup = 'naïve';
[h, C, Cvals] = plotTuningCorrelationsOverReps(v, metric);

v.dataFilter.subjectGroup = 'trained';
[h, C, Cvals] = plotTuningCorrelationsOverReps(v, metric);

v.dataFilter.subjectGroup = 'uncoupled';
[h, C, Cvals] = plotTuningCorrelationsOverReps(v, metric);



% all stimulus dPCs
v.dataFilter.mode_name = 'dpca';
v.dataFilter.mode_OI = 'all_stimulus';
v.dataFilter.mode_method = 'isolate';

v.dataFilter.subjectGroup = 'naïve';
v.dataFilter.mode_file = 'dpca_naive.mat';
[h, C, Cvals] = plotTuningCorrelationsOverReps(v, metric);

v.dataFilter.subjectGroup = 'trained';
v.dataFilter.mode_file = 'dpca_trained.mat';
[h, C, Cvals] = plotTuningCorrelationsOverReps(v, metric);

v.dataFilter.subjectGroup = 'uncoupled';
v.dataFilter.mode_file = 'dpca_uncoupled.mat';
[h, C, Cvals] = plotTuningCorrelationsOverReps(v, metric);



% stimulus dPCs except #1
v.dataFilter.mode_name = 'dpca';
v.dataFilter.mode_OI = 'stimulus';
v.dataFilter.mode_method = 'isolate';

v.dataFilter.subjectGroup = 'naïve';
v.dataFilter.mode_file = 'dpca_naive.mat';
[h, C, Cvals] = plotTuningCorrelationsOverReps(v, metric);

v.dataFilter.subjectGroup = 'trained';
v.dataFilter.mode_file = 'dpca_trained.mat';
[h, C, Cvals] = plotTuningCorrelationsOverReps(v, metric);

v.dataFilter.subjectGroup = 'uncoupled';
v.dataFilter.mode_file = 'dpca_uncoupled.mat';
[h, C, Cvals] = plotTuningCorrelationsOverReps(v, metric);

