%% Knobs

s = false; % save figures to files?

% for sliding windows
window_duration = 1; % [seconds]
t_lim_sec = [-5 35]; % from 5 sec before to 35 seconds after stimulus onset
overlap = .3; % [seconds]

%% Load dataset
filename = 'odorexp004_IC1_130625.mat';
experiment = load(filename).a; % Experiment object

%% initialize output figure saving
cfg = PlotConfig('theme', 'light');

figs = FigureSaver;
figs.outputfolder = fullfiletol('figures',extractBefore(filename,'.'));
figs.outputfile = 'all_plots.pdf';
figs.config = cfg;

%% Plotting average similarity matrices and related metrics for each experimental group.
v = ExperimentViewer(experiment);
v.dataFilter.traceType = 'dFoverF_good';
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
hf = v.plotRepetitionDistances([1 20],'correlation'); % outputs 2 figures
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
hf = v.plotDiscriminationPerformanceMats([1 20], 'correlation','all trials',1:5,false); % outputs 1 figure
figs.title = 'Template-match performance comparison';
figs.append(hf);
close(hf)

% template matching, second by second
window_duration = 1; % [seconds]
t_lim_sec = [-1 40]; % from 1 sec before to 40 seconds after stimulus onset
nwindows = ceil(diff(t_lim_sec)/window_duration)+1;
times = linspace(t_lim_sec(1),t_lim_sec(2),nwindows); % 1 per sec
windows = times(:) + [0 1];
out = cell(nwindows,1);
for i = 1:nwindows
    [hf, out{i}] = v.plotDiscriminationPerformanceMats(windows(i,:), 'correlation','all trials',1:5,false); % outputs 1 figure
    figs.title = ['Template-match: sec', num2str(windows(i,1)), '-', num2str(windows(i,2))];
    figs.append(hf);
    close(hf)
end

% template matching, stimulus window, focus on performance for novel
% stimuli
hf = v.plotDiscriminationPerformanceMats([1 20], 'correlation','all novel',1:5,false); % outputs 1 figure
figs.title = 'Template-match performance on novel stimuli';
figs.append(hf);
close(hf)

% template matching, stimulus window, focus on performance for Leu
hf = v.plotDiscriminationPerformanceMats([1 20], 'correlation',{'Leu'},1:5,false); % outputs 1 figure
figs.title = 'Template-match performance on Leu';
figs.append(hf);
close(hf)

% template matching, stimulus window, focus on performance for familiar
% stimuli
hf = v.plotDiscriminationPerformanceMats([1 20], 'correlation','all familiar',1:5,false); % outputs 1 figure
figs.title = 'Template-match performance on familiar stimuli';
figs.append(hf);
close(hf)

% template matching, stimulus window, focus on performance for CS+
% stimuli
hf = v.plotDiscriminationPerformanceMats([1 20], 'correlation','all CS+',1:5,false); % outputs 1 figure
figs.title = 'Template-match performance on CS+ stimuli';
figs.append(hf);
close(hf)

% template matching, stimulus window, focus on performance for CS-
% stimuli
hf = v.plotDiscriminationPerformanceMats([1 20], 'correlation','all CS-',1:5,false); % outputs 1 figure
figs.title = 'Template-match performance on CS- stimuli';
figs.append(hf);
close(hf)

% template matching based only on trials 1:4, stimulus window
hf = v.plotDiscriminationPerformanceMats([1 20], 'correlation','all trials',1:4,false); % outputs 1 figure
figs.title = 'Template-match performance comparison';
figs.append(hf);
close(hf)

% template matching based only on trials 2:4, stimulus window
hf = v.plotDiscriminationPerformanceMats([1 20], 'correlation','all trials',2:4,false); % outputs 1 figure
figs.title = 'Template-match performance comparison';
figs.append(hf);
close(hf)

% template matching based only on trials 2:5, stimulus window
hf = v.plotDiscriminationPerformanceMats([1 20], 'correlation','all trials',2:5,false); % outputs 1 figure
figs.title = 'Template-match performance comparison';
figs.append(hf);
close(hf)