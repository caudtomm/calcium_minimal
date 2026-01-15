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


%% FIGURE 2

% plot for naive fish
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'naïve';
hf = figure;
v.plotDistancesHead;
v.dataFilter = dft;
figs.title = 'Naive Group';
figs.append(hf);
close(hf)

% native units repetition correlations for same and different stimuli
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'naïve';
hf = figure;
C = v.plotDistancesHead('plotType','repetitions');
idx = repmat(~triu(ones(size(C.distMat3d,1))),1,1,size(C.distMat3d,3));
C_vals_same = 1-C.distMat3d(idx);
idx(:,1,:) = false(size(idx,1),1,size(idx,3));
C_vals_same_noFTE = 1-C.distMat3d(idx);
idx = false(size(idx));
idx(2:end,1,:) = true(size(idx,1)-1,1,size(idx,3));
C_vals_same_FTE = 1-C.distMat3d(idx);
hf = figure;
C = v.plotDistancesHead('plotType','diff_stimulus_repetitions');
C_vals_diff = 1-C.distMat3d(:);
hf = figure;
cdfplot(C_vals_same);hold on
cdfplot(C_vals_same_noFTE)
cdfplot(C_vals_same_FTE)
cdfplot(C_vals_diff)
box off; grid off, axis square
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
set(gcf, 'color', cfg.bgcol); 
legend({'same stimulus','same stimulus (2>5)','same stimulus (1)','different stimuli'})
[~,p] = kstest2(C_vals_same,C_vals_diff);
disp(['2-sample KS test for ''same'' vs ''diff.'' stimuli - pval: ',num2str(p)])
[~,p] = kstest2(C_vals_same_noFTE,C_vals_same_FTE);
disp(['2-sample KS test for ''no FTE'' vs ''FTE'' (same stimulus) - pval: ',num2str(p)])
v.dataFilter = dft;

% example traces for naive fish #6 (Trp and Leu)
% sort by avg intensity of native units on rep 1
v.dataFilter = dft;
sid = v.subjectTab.name(v.subjectTab.group=="naïve"); sid = sid(6);
odors = {'Trp','Leu'};
v.dataFilter.repetitions = [1 3 5];
v.dataFilter.interval = [-5 35];
v.dataFilter.subjectIDs = sid;
for i=1:numel(odors)
    thisodor = odors{i};
    v.dataFilter.stims_allowed = {thisodor};
    
    % native units
    v.dataFilter.mode_name = 'native_units';
    v.dataFilter.mode_method = 'mode_values';
    v.dataFilter.mode_OI = 'all';
    v.dataFilter.mode_file = '';
    [~,out] = v.plotExampleTraces;
    idx = out.idx(1);
    title([thisodor,' - ','units'])

    % isolate stimulus dPCs
    v.dataFilter.mode_name = 'dpca';
    v.dataFilter.mode_OI = 'stimulus';
    v.dataFilter.mode_method = 'isolate';
    v.dataFilter.mode_file = 'dpca_naive.mat';
    v.plotExampleTraces(idx);
    title([thisodor,' - ','isolate stim dPCs'])

    % subtract stimulus dPCs
    v.dataFilter.mode_method = 'subtract';
    v.plotExampleTraces(idx);
    title([thisodor,' - ','subtract stim dPCs'])

end


% dPCA decomposition
v.dataFilter = dft;
v.dataFilter.mode_name = 'dpca';
v.dataFilter.mode_OI = 'all';
v.dataFilter.mode_method = 'mode_values';
v.dataFilter.subjectGroup = 'naïve';
v.dataFilter.mode_file = 'dpca_naive.mat';
hf = figure;
outMat_dn = v.plotDistancesHead;
v.dataFilter.mode_OI = 'all_stimulus';
hf = figure;
outMat_up = v.plotDistancesHead;
v.dataFilter.mode_OI = 'non-stimulus';
hf = figure;
outMat_up = v.plotDistancesHead;
v.dataFilter.mode_OI = 'novelty';
v.dataFilter.interval = [0 3];
hf = figure;
c = v.plotTrialActivityMetricHead('method','avg intensity', 'plotType', 'boxplot_repetitions');
%
v.dataFilter = dft;


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


% stimulus dPC weights
v.dataFilter = dft;
v.dataFilter.mode_name = 'dpca';
v.dataFilter.mode_OI = 'stimulus';
v.dataFilter.mode_method = 'mode_values';
v.dataFilter.subjectGroup = 'naïve';
v.dataFilter.mode_file = 'dpca_naive.mat';
m = ModeSelector(v).extract;
c = []; % [w1 N; w2 N; ...]
for i= 1:numel(m.coeffs)
    thisdata = m.coeffs{i}(:,m.parseModeOI(i));
    thisdata = sum(thisdata);
    c = [c; thisdata(:)];
end
%histogram
figure; histogram(c, 100, 'FaceColor','k','EdgeAlpha',0);
box off; axis square
axis tight
xlabel('sum of weights'); ylabel('histogram')
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
set(gcf, 'color', cfg.bgcol);


% stimulus dPC weights
v.dataFilter = dft;
v.dataFilter.mode_name = 'dpca';
v.dataFilter.mode_OI = 'stimulus';
v.dataFilter.mode_method = 'mode_values';
v.dataFilter.subjectGroup = 'naïve';
v.dataFilter.mode_file = 'dpca_naive.mat';
m = ModeSelector(v).extract;
c = [];
for i= 1:numel(m.coeffs)
    thisdata = m.coeffs{i}(:,m.parseModeOI(i));
    thisdata = sum(thisdata);
    c = [c; thisdata(:)];
end
%histogram
figure; histogram(c, 100, 'FaceColor','k','EdgeAlpha',0);
box off; axis square
axis tight
xlabel('sum of weights'); ylabel('histogram')
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
set(gcf, 'color', cfg.bgcol);
% 
c = []; % [w1 N; w2 N; ...]
for i= 1:numel(m.coeffs)
    thisdata = m.coeffs{i}(:,m.parseModeOI(i));
    c = [c; thisdata(:), ones(numel(thisdata),1)*v.filtered_traces{i}.N];
end
%histogram
figure; histogram(abs(c(:,1)).*c(:,2), 100, 'FaceColor','k','EdgeAlpha',0);
box off; axis square
yscale log
xscale log
axis tight
xlabel('|w|*N'); ylabel('histogram')
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
set(gcf, 'color', cfg.bgcol);
% shower plot
figure; scatter(c(:,2),abs(c(:,1)),5,'k','filled')
xlim([50 350])
hold on
x = min(xlim):1:max(xlim);
y = 1./x;
plot(x,y,'--r','LineWidth',1)
yscale log
for i= 1:numel(m.coeffs)
    avg_val = mean(abs(m.coeffs{i}(:,m.parseModeOI(i))),"all","omitmissing");
    N = v.filtered_traces{i}.N;
    scatter(N,avg_val,35,'filled','Marker','diamond','MarkerFaceColor','r')
end
xlabel('N')
ylabel('|w|')
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
set(gcf, 'color', cfg.bgcol);
% mode value for each stimulus in 5 example fish
nsubjects = numel(v.filtered_traces);
idxs = randperm(nsubjects,5);
hf= figure;
n=1;
for i=1:5
    idx = idxs(i);
    thismodevals = m.values{idx};
    [N,nmodes,ntrials] = size(thismodevals);
    avg_vals = squeeze(mean(reshape(thismodevals,[N,nmodes,ntrials/6,6]),[1,3],'omitmissing')); % avg value for each stimulus
    Y = pdist(avg_vals,'cosine'); % figure; imagesc(squareform(Y));
    Z = linkage(Y,'average');
    cophenet(Z,Y)
    T = cluster(Z,"maxclust",6);
    [T, idx] = sort(T);
    subplot(5,2,n); n = n+1; imagesc(zscore([avg_vals(idx,:)]'))
    axis image
    xticks(1:nmodes); xticklabels(T)
    ylabel('stimuli'); xlabel('dPC cluster')
    colormap(cfg.colormapName)
    clim([-2 2])
    set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
    subplot(5,2,n); n = n+1; dendrogram(Z,'Reorder',idx)
    set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
end
set(gcf, 'color', cfg.bgcol);
% for each odor, how many dPCs at least 1 STD above their own baseline?
allportion = [];
th_bins = -2:.1:2;
for i_bin = 1:numel(th_bins)
    thisth = th_bins(i_bin);
    thisportion = [];
    for i = 1:nsubjects
        thismodevals = m.values{i};
        [N,nmodes,ntrials] = size(thismodevals);
        avg_vals = squeeze(mean(reshape(thismodevals,[N,nmodes,ntrials/6,6]),[1,3],'omitmissing')); % avg value for each stimulus
        avg_vals = nanzscore(avg_vals,[],2);
        thisportion = [thisportion, sum(avg_vals>=thisth)./nmodes];
    end
    allportion = [allportion;thisportion];
end
allportion = allportion*100; % (%)
figure; 
t = th_bins;
mu = mean(allportion,2);
err = std(allportion,[],2);
% Shaded error bars
fill([t fliplr(t)], [mu - err; flipud(mu + err)]', ...
     'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % grey shade
hold on
plot(t,mu,'k-','LineWidth',1)
box off; axis square tight
xlabel('threshold [z-score]'); ylabel('% dPCs above threshold')
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
set(gcf, 'color', cfg.bgcol);
% correlation of loadings across dPCs (orthogonal?)
c = [];
for i= 1:nsubjects
    thisdata = m.coeffs{i}(:,m.parseModeOI(i));
    thisdata = [1-pdist(thisdata','cosine')]';
    c = [c; thisdata];
end
figure; histogram(c,100,'FaceColor','k','EdgeAlpha',0);
box off; axis square
xlim([-1 1])
axis tight
xlabel('w similarity'); ylabel('histogram')
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
set(gcf, 'color', cfg.bgcol);

% variance explained by each PC/dPC
v.dataFilter = dft;
v.dataFilter.mode_name = 'dpca';
v.dataFilter.mode_OI = 'all_stimulus';
v.dataFilter.mode_method = 'mode_values';
v.dataFilter.subjectGroup = 'naïve';
v.dataFilter.mode_file = 'dpca_naive.mat';
m = ModeSelector(v).extract;
nsubjects = numel(m.fullout);
data = nan(3,20,nsubjects);
for i=1:nsubjects
    thisdata = m.fullout{i}.explVar;
    data(1,:,i) = thisdata.cumulativePCA;
    
    idx = m.parseModeOI(i);
    idx_novelty = find(idx,1);
    idx(idx_novelty)=false;

    data(2,1:sum(idx),i) = cumsum(thisdata.componentVar(idx));
    data(3,1,i) = thisdata.componentVar(idx_novelty);
    
end
hf = figure; clear b
t = 1:20;
mu = squeeze(mean(data(1,:,:),3,'omitmissing'))'; % PCs
err = std(data(1,:,:),[],3,'omitmissing')';
fill([t fliplr(t)], [mu - err; flipud(mu + err)]', ...
     'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
hold on
b(1) = plot(t,mu,'k-','LineWidth',1);
mu = squeeze(mean(data(2,:,:),3,'omitmissing'))'; % stimulus dPCs
idx = find(isnan(mu),1);
t = t(1:idx-1);
mu = mu(1:idx-1);
err = std(data(2,1:idx-1,:),[],3,'omitmissing')';
fill([t fliplr(t)], [mu - err; flipud(mu + err)]', ...
     'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
b(2) = plot(t,mu,'r-','LineWidth',1);
mu = squeeze(mean(data(3,1,:),3,'omitmissing'))'; % novelty dPC
err = std(data(3,1,:),[],3,'omitmissing')';
b(3) = scatter(1,mu,50,'blue','filled');
errorbar(mu,err,'Color','blue')
legend(b,{'PCs','dPCs','novelty dPC'})
box off; axis square
xlabel('component #'); ylabel('histogram')
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
set(gcf, 'color', cfg.bgcol);


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


% only novelty dPC
%
% mode weights
v.dataFilter = dft;
v.dataFilter.mode_name = 'dpca';
v.dataFilter.mode_OI = 'novelty';
v.dataFilter.mode_method = 'mode_values';
v.dataFilter.subjectGroup = 'naïve';
v.dataFilter.mode_file = 'dpca_naive.mat';
m = ModeSelector(v).extract;
c = []; % [w1 N; w2 N; ...]
for i= 1:numel(m.coeffs)
    c = [c; m.coeffs{i}(:,m.parseModeOI(i)), ones(v.filtered_traces{i}.N,1)*v.filtered_traces{i}.N];
end
figure; scatter(c(:,2),abs(c(:,1)),5,'filled')
xlim([50 350])
hold on
x = min(xlim):1:max(xlim);
y = 1./x;
plot(x,y,'--r','LineWidth',1)
yscale log
for i= 1:numel(m.coeffs)
    avg_val = mean(abs(m.coeffs{i}(:,m.parseModeOI(i))),"all","omitmissing");
    N = v.filtered_traces{i}.N;
    scatter(N,avg_val,35,'filled','Marker','diamond','MarkerFaceColor','r')
end
xlabel('N')
ylabel('|w|')
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
set(gcf, 'color', cfg.bgcol);
% correlation with general suppression
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'naïve';
supp = v.plotUnitActivityMetricHead('method','general suppression score');
supp = cell2mat(supp);
figure; plotHeatmapAndIsoclines(supp,log(abs(c(:,1))),25,1,0,0)
xlabel('suppression [a.u]'); ylabel('log(|w|)')
axis square
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
set(gcf, 'color', cfg.bgcol); 
figure; scatter(c(:,1),supp,5,'k','filled')
ylabel('suppression [a.u]'); xlabel('w')
ylim([-2.5 2.5])
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
set(gcf, 'color', cfg.bgcol); 
[pc,pval] = corrcoef(c(:,1),supp);
pc = pc(2); pval = pval(2);
coeff_det = pc^2;
disp(['Pearson Correlation: ', num2str(pc)])
disp(['p value: ', num2str(pval)])
disp(['Coeff. of Determination (variance explained): ', num2str(coeff_det*100), '%'])
disp('')
p = polyfit(c(:,1),supp, 1); % trend line
x = mean(c(:,1),'omitmissing') + std(c(:,1),'omitmissing') .* [-1 1];
hold on; plot(x,polyval(p,x),'r','LineWidth',2)
% correlation with odor-specific suppression
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'naïve';
supp = v.plotUnitActivityMetricHead('method','stimulus specific suppression score');
supp = cell2mat(supp);
figure; plotHeatmapAndIsoclines(supp,log(abs(c(:,1))),25,1,0,0)
xlabel('suppression [a.u]'); ylabel('log(|w|)')
axis square
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
set(gcf, 'color', cfg.bgcol); 
figure; scatter(c(:,1),supp,5,'k','filled')
ylabel('suppression [a.u]'); xlabel('w')
ylim([-2.5 2.5])
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
set(gcf, 'color', cfg.bgcol); 
[pc,pval] = corrcoef(c(:,1),supp);
pc = pc(2); pval = pval(2);
coeff_det = pc^2;
disp(['Pearson Correlation: ', num2str(pc)])
disp(['p value: ', num2str(pval)])
disp(['Coeff. of Determination (variance explained): ', num2str(coeff_det*100), '%'])
disp('')
p = polyfit(c(:,1),supp, 1); % trend line
x = mean(c(:,1),'omitmissing') + std(c(:,1),'omitmissing') .* [-1 1];
hold on; plot(x,polyval(p,x),'r','LineWidth',2)
% correlation with population intensity
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'naïve';
v.dataFilter.repetitions = 1;
unitint = v.plotUnitActivityMetricHead('method','avg intensity');
unitint = mean(cell2mat(unitint), 2, 'omitmissing');
figure; plotHeatmapAndIsoclines(unitint,log(abs(c(:,1))),25,1,0,0)
xlabel('iFR on rep 1 [Hz]'); ylabel('log(|w|)')
axis square
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
set(gcf, 'color', cfg.bgcol); 
figure; scatter(c(:,1),unitint,5,'k','filled')
ylabel('iFR on rep 1 [Hz]'); xlabel('w')
ylim([0 .8])
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
set(gcf, 'color', cfg.bgcol, 'Position', [100 100 433 104]); 
[pc,pval] = corrcoef(c(:,1),unitint);
pc = pc(2); pval = pval(2);
coeff_det = pc^2;
disp(['Pearson Correlation: ', num2str(pc)])
disp(['p value: ', num2str(pval)])
disp(['Coeff. of Determination (variance explained): ', num2str(coeff_det*100), '%'])
disp('')
p = polyfit(c(:,1),unitint, 1); % trend line
x = mean(c(:,1),'omitmissing') + std(c(:,1),'omitmissing') .* [-1 1];
hold on; plot(x,polyval(p,x),'r','LineWidth',2)
% histogram
figure; histogram(c(:,1),100, 'FaceColor','k','EdgeAlpha',0);
box off; axis square
xlabel('w'); ylabel('histogram')
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
set(gcf, 'color', cfg.bgcol);

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



%% FIGURE 3


% stimulus dPCs except #1 (trained)
v.dataFilter = dft;
v.dataFilter.mode_name = 'dpca';
v.dataFilter.mode_OI = 'stimulus';
v.dataFilter.mode_method = 'subtract';
v.dataFilter.subjectGroup = 'trained';
v.dataFilter.mode_file = 'dpca_trained.mat';
hf = figure;
outMat_dn = v.plotDistancesHead;
v.dataFilter.mode_method = 'isolate';
hf = figure;
outMat_up = v.plotDistancesHead;
v.dataFilter.mode_method = 'mode_values';
hf = figure;
outMat_up = v.plotDistancesHead;
%
v.dataFilter = dft;


% stimulus dPCs except #1 (uncoupled)
v.dataFilter = dft;
v.dataFilter.mode_name = 'dpca';
v.dataFilter.mode_OI = 'stimulus';
v.dataFilter.mode_method = 'subtract';
v.dataFilter.subjectGroup = 'uncoupled';
v.dataFilter.mode_file = 'dpca_uncoupled.mat';
hf = figure;
outMat_dn = v.plotDistancesHead;
v.dataFilter.mode_method = 'isolate';
hf = figure;
outMat_up = v.plotDistancesHead;
v.dataFilter.mode_method = 'mode_values';
hf = figure;
outMat_up = v.plotDistancesHead;
%
v.dataFilter = dft;

% plot interrep correlations and delta for trained fish
v.dataFilter = dft;
hf = figure;
v.dataFilter.subjectGroup = 'trained';
v.plotDistancesHead();
hf = figure;
v.dataFilter.subjectGroup = 'naïve';
naive3d = v.plotDistancesHead('plotType','repetitions');
v.dataFilter.subjectGroup = 'trained';
trained3d = v.plotDistancesHead('plotType','repetitions');
naive2d = 1-mean(naive3d.distMat3d,3,'omitmissing');
trained2d = 1-mean(trained3d.distMat3d,3,'omitmissing');
% lij
inMat = trained2d-naive2d;
triu_idx = logical(triu(true(size(inMat))));
h = hf.Children(2).Children;
set(h,'AlphaData',[~triu_idx]')
hf = figure;
h = imagesc(inMat,'AlphaData',~triu_idx);
clim([-.2 .2])
colormap('jet')
n_repetitions = height(h.CData);
axis square; hold on
xticks(1:n_repetitions); xticklabels(1:n_repetitions); xtickangle(90)
yticks(1:n_repetitions); yticklabels(1:n_repetitions)
xlabel('Repetition')
ylabel('Repetition')
colormap(cfg.colormapName)
a = colorbar('Color',cfg.axcol);
a.Label.String = 'delta correlation';
a.Label.FontSize= gca().FontSize;
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
set(gcf, 'color', cfg.bgcol); 
hold off

% plot interrep correlations and delta for uncoupled fish
v.dataFilter = dft;
hf = figure;
v.dataFilter.subjectGroup = 'uncoupled';
v.plotDistancesHead();
hf = figure;
v.dataFilter.subjectGroup = 'naïve';
naive3d = v.plotDistancesHead('plotType','repetitions');
v.dataFilter.subjectGroup = 'uncoupled';
trained3d = v.plotDistancesHead('plotType','repetitions');
naive2d = 1-mean(naive3d.distMat3d,3,'omitmissing');
trained2d = 1-mean(trained3d.distMat3d,3,'omitmissing');
% lij
inMat = trained2d-naive2d;
triu_idx = logical(triu(true(size(inMat))));
h = hf.Children(2).Children;
set(h,'AlphaData',[~triu_idx]')
hf = figure;
h = imagesc(inMat,'AlphaData',~triu_idx);
clim([-.2 .2])
colormap('jet')
n_repetitions = height(h.CData);
axis square; hold on
xticks(1:n_repetitions); xticklabels(1:n_repetitions); xtickangle(90)
yticks(1:n_repetitions); yticklabels(1:n_repetitions)
xlabel('Repetition')
ylabel('Repetition')
colormap(cfg.colormapName)
a = colorbar('Color',cfg.axcol);
a.Label.String = 'delta correlation';
a.Label.FontSize= gca().FontSize;
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
set(gcf, 'color', cfg.bgcol); 
hold off


% familiar (dn) vs novel (up) odors (trained fish)
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'trained';
v.dataFilter.stims_allowed = 'all novel';
hf = figure;
v.plotDistancesHead('plotType','repetitions');
h = hf.Children(2).Children;
n_repetitions = height(h.CData);
triu_idx = logical(triu(true(n_repetitions)));
set(h,'AlphaData',[~triu_idx]');
v.dataFilter.stims_allowed = 'all familiar';
hf = figure;
v.plotDistancesHead('plotType','repetitions');
h = hf.Children(2).Children;
set(h,'AlphaData',~triu_idx);

% familiar (dn) vs novel (up) odors (uncoupled fish)
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'uncoupled';
v.dataFilter.stims_allowed = 'all novel';
hf = figure;
v.plotDistancesHead('plotType','repetitions');
h = hf.Children(2).Children;
n_repetitions = height(h.CData);
triu_idx = logical(triu(true(n_repetitions)));
set(h,'AlphaData',[~triu_idx]');
v.dataFilter.stims_allowed = 'all familiar';
hf = figure;
v.plotDistancesHead('plotType','repetitions');
h = hf.Children(2).Children;
set(h,'AlphaData',~triu_idx);
