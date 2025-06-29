%% Knobs

s = true; % save figures to files?

%% Load dataset
filename = 'odorexp004_IC1_130625.mat';
experiment = load(filename).a; % Experiment object

%% initialize output figure saving
cfg = PlotConfig('theme', 'dark');

figs = FigureSaver;
figs.outputfolder = fullfiletol('figures',extractBefore(filename,'.'));
figs.outputfile = 'all_plots.pdf';
figs.config = cfg;

%% Plotting average similarity matrices and related metrics for each experimental group.
v = ExperimentViewer(experiment);
v.plotConfig = cfg;

% plot for naive fish
v.dataFilter.subjectGroup = 'naïve';
hf = v.plotDistances('ps_lim',[1,20],'trial_sorting','stim_id');
figs.title = 'Naive Group';
figs.append(hf);
close(hf)

% plot for trained fish
v.dataFilter.subjectGroup = 'trained';
hf = v.plotDistances('ps_lim',[1,20],'trial_sorting','stim_id');
figs.title = 'Trained Groups';
figs.append(hf);
close(hf)

% plot for uncoupled fish
v.dataFilter.subjectGroup = 'trained';
hf = v.plotDistances('ps_lim',[1,20],'trial_sorting','stim_id');
figs.title = 'Trained Groups';
figs.append(hf);
close(hf)
