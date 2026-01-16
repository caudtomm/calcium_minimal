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
 

%% FIGURE 4

% naive intertrial distances
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'naïve';
hf = figure;
v.plotDistancesHead;
v.dataFilter = dft;

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
