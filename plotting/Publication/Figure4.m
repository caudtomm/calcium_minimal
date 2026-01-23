dbstop if error

s = true; % save figures to files?
savepath = 'bin4';
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
 

%% FIGURE 4


%% intertrial correlations, without pre-stimulus correlations
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'trained';
v.dataFilter.stims_allowed = 'all stimuli';
v.dataFilter.interval = [-22 -2];
hf = figure;
subplot(221)
C = v.plotDistancesHead;
C= 1-C.distMat3d;
hold on
v.dataFilter.trial_sorting = 'chronological';
v.plotDistancesHead;
v.dataFilter.trial_sorting = 'stim_id';
subplot(222); v.dataFilter.interval = [1 20];
C1 = v.plotDistancesHead;
C1 = 1- C1.distMat3d;
hold on
imagesc(mean(C1-C,3,'omitmissing'))
clim([-.2 .2])
subplot(234)
v.dataFilter.interval = [-22 -2];
C = v.plotDistancesHead('plotType','repetitions');
C= 1-C.distMat3d;
subplot(235); v.dataFilter.interval = [1 20];
v.plotDistancesHead('plotType','repetitions');
subplot(236); v.dataFilter.interval = [1 20];
C1 = v.plotDistancesHead('plotType','repetitions');
C1 = 1- C1.distMat3d;
hold on
imagesc(mean(C1-C,3,'omitmissing'))
clim([-.2 .2])
v.dataFilter = dft;


%% intertrial correlations, without pre-stimulus correlations
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'uncoupled';
v.dataFilter.stims_allowed = 'all stimuli';
v.dataFilter.interval = [-22 -2];
[~,base_events] = ModeSelector(v).extract;
v.dataFilter.interval = [1 20];
[~,odor_events,all_labs] = ModeSelector(v).extract;
base_events = cellfun(@(x) mean(x,1,'omitmissing'),base_events,'UniformOutput',false);
odor_events = cellfun(@(x) mean(x,1,'omitmissing'),odor_events,'UniformOutput',false);
diff_events = cellfun(@(x,y) y-x,base_events,odor_events,'UniformOutput',false);
figure; plotDistances(diff_events,'full','correlation',all_labs{1},cfg);
clim([0 .5])
cfg.figSize = "medium";
cfg.setFigure;
cfg.saveFigure(gcf,'uncoupled stim-prestim intertrial corr', saveType)
v.dataFilter = dft;

figure; plotDistances(diff_events,'repetitions','correlation',all_labs{1},cfg);
clim([0 .5])
cfg.figSize = "small";
cfg.setFigure;
cfg.saveFigure(gcf,'uncoupled stim-prestim intertrial corr repetitions', saveType)
v.dataFilter = dft;


%% intertrial distances (novel or familiar stimuli)
%
stimclass = 'novel';

% naive
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'naïve';
plotIntertrialDist(v,stimclass,'naive',cfg,saveType);
% trained
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'trained';
plotIntertrialDist(v,stimclass,'trained',cfg,saveType);
% uncoupled
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'uncoupled';
plotIntertrialDist(v,stimclass,'uncoupled',cfg,saveType);

function plotIntertrialDist(v,stimclass,grouptag,cfg,saveType)
    if strcmp(stimclass,'familiar')
        v.dataFilter.stims_allowed = {'Arg','Ala','His'};
    elseif strcmp(stimclass,'novel')
        v.dataFilter.stims_allowed = {'Trp','Ser','Leu'};
    end
    hf = figure;
    v.plotDistancesHead;
    xticks([]); yticks([]); xlabel(''); ylabel(''); title('')
    cfg.setFigure;
    cfg.saveFigure(gcf,[grouptag,' intertrial corr ',stimclass], saveType)
end

%% influence of valence on intertrial correlations

% trained1: CS+ (Arg + Ala)
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'trained1';
v.dataFilter.stims_allowed = 'all CS+';
hf = figure;
subplot(311);
C1 = v.plotDistancesHead('plotType','diff_stimulus_repetitions');
v.dataFilter = dft;
xticks([]); yticks([]); xlabel(''); ylabel(''); title('trained1')
% trained2: same (Arg + Ala) <- CS- vs CS+
v.dataFilter = dft;
subplot(312);
v.dataFilter.subjectGroup = 'trained2';
v.dataFilter.stims_allowed = {'Arg','Ala'};
C2 = v.plotDistancesHead('plotType','diff_stimulus_repetitions');
v.dataFilter = dft;
xticks([]); yticks([]); xlabel(''); ylabel(''); title('trained2')
subplot(313);
C1 = 1 - C1.distMat3d(:);
C2 = 1 - C2.distMat3d(:);
y = [[C1 ones(length(C1),1)] ; [C2 2*ones(length(C2),1)]];
boxplot(y(:,1),y(:,2)); axis square
hold on; scatter([1,2],[mean(C1,'omitmissing'),mean(C2,'omitmissing')],15,'r','filled')
xticks([1 2]);xticklabels({'trained1','trained2'}); ylabel('r')
ylim([-.2 1])
[~,p] = ttest2(C1,C2);
cfg.figSize = "small";
cfg.aspRatioType = 'tall';
cfg.setFigure;
cfg.saveFigure(gcf,'trained1 vs trained2 Arg-Ala intertrial corr', saveType)
cfg = v.plotConfig;

% trained1: CS+ (Ala + His)
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'trained1';
v.dataFilter.stims_allowed = {'Ala','His'};
hf = figure;
subplot(311);
C1 = v.plotDistancesHead('plotType','diff_stimulus_repetitions');
v.dataFilter = dft;
xticks([]); yticks([]); xlabel(''); ylabel(''); title('trained1')
% trained2: same (Arg + Ala) <- CS- vs CS+
v.dataFilter = dft;
subplot(312);
v.dataFilter.subjectGroup = 'trained2';
v.dataFilter.stims_allowed = 'all CS+';
C2 = v.plotDistancesHead('plotType','diff_stimulus_repetitions');
v.dataFilter = dft;
xticks([]); yticks([]); xlabel(''); ylabel(''); title('trained2')
subplot(313);
C1 = 1 - C1.distMat3d(:);
C2 = 1 - C2.distMat3d(:);
y = [[C1 ones(length(C1),1)] ; [C2 2*ones(length(C2),1)]];
boxplot(y(:,1),y(:,2)); axis square
hold on; scatter([1,2],[mean(C1,'omitmissing'),mean(C2,'omitmissing')],15,'r','filled')
xticks([1 2]);xticklabels({'trained1','trained2'}); ylabel('r')
ylim([-.2 1])
[~,p] = ttest2(C1,C2);
cfg.figSize = "small";
cfg.aspRatioType = 'tall';
cfg.setFigure;
cfg.saveFigure(gcf,'trained1 vs trained2 Ala-His intertrial corr', saveType)
cfg = v.plotConfig;

% CS+ vsd CS+
v.dataFilter = dft;
hf = figure;
subplot(311);
v.dataFilter.subjectGroup = 'trained1';
v.dataFilter.stims_allowed = {'Arg','Ala'};
C1_1 = v.plotDistancesHead('plotType','diff_stimulus_repetitions');
v.dataFilter.subjectGroup = 'trained2';
v.dataFilter.stims_allowed = {'Ala','His'};
C1_2 = v.plotDistancesHead('plotType','diff_stimulus_repetitions');
n1 = size(C1_1.distMat3d,3);
n2 = size(C1_2.distMat3d,3);
C1 = zeros(5,5,n1+n2);
C1(:,:,1:n1) = 1- C1_1.distMat3d;
C1(:,:,(n1+1):end) = 1- C1_2.distMat3d;
hold on
imagesc(mean(C1,3,'omitmissing'));
xticks([]); yticks([]); xlabel(''); ylabel(''); title('CS+ vs CS+')
% CS+ vs CS-
v.dataFilter = dft;
subplot(312);
v.dataFilter.subjectGroup = 'trained1';
v.dataFilter.stims_allowed = {'Ala','His'};
C2_1 = v.plotDistancesHead('plotType','diff_stimulus_repetitions');
v.dataFilter.subjectGroup = 'trained2';
v.dataFilter.stims_allowed = {'Arg','Ala'};
C2_2 = v.plotDistancesHead('plotType','diff_stimulus_repetitions');
n1 = size(C2_1.distMat3d,3);
n2 = size(C2_2.distMat3d,3);
C2 = zeros(5,5,n1+n2);
C2(:,:,1:n1) = 1- C2_1.distMat3d;
C2(:,:,(n1+1):end) = 1- C2_2.distMat3d;
hold on
imagesc(mean(C2,3,'omitmissing'));
xticks([]); yticks([]); xlabel(''); ylabel(''); title('CS+ vs CS-')
subplot(313);
C1 = C1(:); C2 = C2(:);
y = [[C1 ones(length(C1),1)] ; [C2 2*ones(length(C2),1)]];
boxplot(y(:,1),y(:,2)); axis square
hold on; scatter([1,2],[mean(C1,'omitmissing'),mean(C2,'omitmissing')],5,'r','filled')
xticks([1 2]);xticklabels({'CS+ vs CS+','CS+ vs CS-'}); ylabel('r')
ylim([-.2 1])
[~,p] = ttest2(C1,C2);
cfg.figSize = "small";
cfg.aspRatioType = 'tall';
cfg.setFigure;
cfg.saveFigure(gcf,'same vs diff valence intertrial corr', saveType)
cfg = v.plotConfig;


%%
v.dataFilter = dft;
stimclass = 'novel';
% v.dataFilter.stims_allowed = 'all stimuli';
% v.dataFilter.stims_allowed = {'Arg','Ala','His'};
v.dataFilter.stims_allowed = {'Trp','Ser','Leu'};


% native units repetition correlations for same and different stimuli
v.dataFilter.subjectGroup = 'naïve';
C_naive = prova(v,saveType,'naive','',stimclass);
v.dataFilter.subjectGroup = 'trained';
C_trained = prova(v,saveType,'trained','',stimclass);
v.dataFilter.subjectGroup = 'uncoupled';
C_uncoupled = prova(v,saveType,'uncoupled','',stimclass);
close all

y = [C_naive;C_trained;C_uncoupled];
g = [repmat({'naive'},numel(C_naive),1); ...
    repmat({'trained'},numel(C_trained),1); ...
    repmat({'uncoupled'},numel(C_uncoupled),1)];
hf = figure;
boxplot(y,g);
ylim([-1 1]); ylabel('r')
[p,~,stats] = kruskalwallis(y,g,'off');
figure;multcompare(stats,'Display','on')

%
% native units repetition correlations for same and different stimuli (early peak interval)
v.dataFilter.interval = [.5 3];
v.dataFilter.subjectGroup = 'naïve';
C_naive = prova(v,saveType,'naive','initpeak',stimclass);
v.dataFilter.subjectGroup = 'trained';
C_trained = prova(v,saveType,'trained','initpeak',stimclass);
v.dataFilter.subjectGroup = 'uncoupled';
C_uncoupled = prova(v,saveType,'uncoupled','initpeak',stimclass);
close all

%
% native units repetition correlations for same and different stimuli (main response body interval)
v.dataFilter.interval = [4 20];
v.dataFilter.subjectGroup = 'naïve';
C_naive = prova(v,saveType,'naive','body',stimclass);
v.dataFilter.subjectGroup = 'trained';
C_trained = prova(v,saveType,'trained','body',stimclass);
v.dataFilter.subjectGroup = 'uncoupled';
C_uncoupled = prova(v,saveType,'uncoupled','body',stimclass);
close all

%
% native units repetition correlations for same and different stimuli (early peak interval)
v.dataFilter.interval = [20.5 25];
v.dataFilter.subjectGroup = 'naïve';
prova(v,saveType,'naive','off',stimclass)
v.dataFilter.subjectGroup = 'trained';
prova(v,saveType,'trained','off',stimclass)
v.dataFilter.subjectGroup = 'uncoupled';
prova(v,saveType,'uncoupled','off',stimclass)


function C_vals_same = prova(v, saveType, groupname, tag,stimclass)
    if nargin<4; tag = ''; end
    cfg = v.plotConfig;
    hf = figure;
    subplot(311);
    C = v.plotDistancesHead('plotType','repetitions');
    xticks([]); yticks([]); xlabel(''); ylabel(''); title('')
    idx = repmat(~triu(ones(size(C.distMat3d,1))),1,1,size(C.distMat3d,3));
    C_vals_same = 1-C.distMat3d(idx);
    idx(:,1,:) = false(size(idx,1),1,size(idx,3));
    C_vals_same_noFTE = 1-C.distMat3d(idx);
    idx = false(size(idx));
    idx(2:end,1,:) = true(size(idx,1)-1,1,size(idx,3));
    C_vals_same_FTE = 1-C.distMat3d(idx);
    idx = logical(repmat(triu(ones(size(C.distMat3d,1)),1)-triu(ones(size(C.distMat3d,1)),2),1,1,size(C.distMat3d,3)));
    C_vals_adj = 1-C.distMat3d(idx);
    C_vals_adj = reshape(C_vals_adj,size(C.distMat3d,1)-1,[])';
    subplot(312);
    thisgroup = v.dataFilter.subjectGroup;
    v.dataFilter.subjectGroup = 'naïve';
    C_naive = v.plotDistancesHead('plotType','repetitions');
    C = mean(C_naive.distMat3d,3,'omitmissing') - mean(C.distMat3d,3,'omitmissing');
    hold on; imagesc(C);
    clim([-.2 .2])
    subplot(313);
    C = v.plotDistancesHead('plotType','diff_stimulus_repetitions');
    xticks([]); yticks([]); xlabel(''); ylabel(''); title('')
    C_vals_diff = 1-C.distMat3d(:);
    cfg.setFigure;
    cfg.saveFigure(gcf,[groupname, ' intertrial corr ',stimclass,' pooled',tag], saveType)
    cfg = v.plotConfig;

    hf = figure;
    mu = mean(C_vals_adj,1,'omitmissing'); err = std(C_vals_adj,[],1,'omitmissing');
    plotLineNShade(1:4, mu,err,'k', cfg)
    xticks([1:4]); ylim([0 1]); xlabel('Repetition #'); ylabel('r')
    cfg.figSize = 'tiny';
    cfg.aspRatioType = 'square';
    cfg.setFigure;
    cfg.saveFigure(gcf,[groupname, ' intertrial corr ',stimclass,' pooled adjecent',tag], saveType)
    cfg = v.plotConfig;

    hf = figure;
    cdfplot(C_vals_same);hold on
    cdfplot(C_vals_same_noFTE)
    cdfplot(C_vals_same_FTE)
    cdfplot(C_vals_diff)
    legend({'same stimulus','same stimulus (2>5)','same stimulus (1)','different stimuli'})
    xlabel('r'); ylabel('CDF'); title('')
    xlim([-.2 1])
    view([90 -90])
    cfg.figSize = 'tiny';
    cfg.aspRatioType = 'tall';
    cfg.setFigure;
    cfg.saveFigure(gcf,[groupname, ' intertrial corr ',stimclass,' pooled CDF',tag], saveType)
    cfg = v.plotConfig;

    [~,p] = kstest2(C_vals_same,C_vals_diff);
    disp(['2-sample KS test for ''same'' vs ''diff.'' stimuli - pval: ',num2str(p)])
    [~,p] = kstest2(C_vals_same_noFTE,C_vals_same_FTE);
    disp(['2-sample KS test for ''no FTE'' vs ''FTE'' (same stimulus) - pval: ',num2str(p)])
    
end

%%
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
