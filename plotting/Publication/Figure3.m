dbstop if error

s = true; % save figures to files?
savepath = 'bin3';
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
cfg.custom.crange = [.1 .7];
cfg.savePath = savepath;

v = ExperimentViewer(experiment);
v.plotConfig = cfg;
cfg = v.plotConfig;

v.dataFilter.traceType = 'pSpike';
v.dataFilter.interval = [1,20];
v.dataFilter.trial_sorting = 'stim_id';
dft = v.dataFilter;
 

%% FIGURE 3


%% unit firing distributions

baseline_interval = [-22 -2];
odor_interval = dft.interval;
yrange = [-.2 .3];
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'trained';
v.dataFilter.trial_sorting = 'chronological';
% over trials
v.dataFilter.stims_allowed = 'all trials';
out = v.plotUnitActivityMetricHead('method','avg intensity');
v.dataFilter.interval = baseline_interval;
outbase = v.plotUnitActivityMetricHead('method','avg intensity');
data = cell2mat(out) - cell2mat(outbase);
[~,~, labs] = ModeSelector(v).extract; labs = labs{1};
hf = figure; b = prettyBoxplot(data,labs,'scatterSize',5,'plotLine',true);
ylabel('Cellwise delta iFR (Hz)')
ylim(yrange)
cfg.figSize = "large";
cfg.aspRatioType = "wide";
cfg.lineWidth = .5;
cfg.setFigure
cfg.saveFigure(gcf,'trained1 unit delta firing distribution', saveType)
cfg = v.plotConfig;
% overall
v.dataFilter.stims_allowed = 'all stimuli';
v.dataFilter.subjectGroup = 'trained1';
v.dataFilter.interval = odor_interval;
out = v.plotUnitActivityMetricHead('method','avg intensity');
v.dataFilter.interval = baseline_interval;
outbase = v.plotUnitActivityMetricHead('method','avg intensity');
data_naive = cell2mat(out) - cell2mat(outbase);
v.dataFilter.subjectGroup = 'trained';
v.dataFilter.interval = odor_interval;
out = v.plotUnitActivityMetricHead('method','avg intensity');
v.dataFilter.interval = baseline_interval;
outbase = v.plotUnitActivityMetricHead('method','avg intensity');
data_trained = cell2mat(out) - cell2mat(outbase);
[~,~,labs] = ModeSelector(v).extract;
labs = labs{1};
repnum = arrayfun(@(i) sum(strcmp(labs(1:i), labs(i))), 1:numel(labs))';
hf = figure;
[stims,~,stimidx] = unique(labs,'stable');
cols = colorcube(10);
hold on;
for i = 1:numel(stims)
    idx = stimidx == i;
    scatter(mean(data_naive(:,idx),1,'omitmissing'), mean(data_trained(:,idx),1,'omitmissing'),...
        30, 'filled', 'CData', cols(i,:), 'MarkerFaceAlpha', 'flat', ...
        'AlphaData', .7 + .3./repnum(idx), 'DisplayName', stims{i});
end
maxrange = max([max(xlim),max(ylim)]);
xlim([0 maxrange]); ylim([0 maxrange]);
plot([0 maxrange],[0 maxrange],'r-','DisplayName','x=y')
u = legendUnq; legend(u)
hold off;
axis square
xlabel('naive iFR (Hz)');
ylabel('trained iFR (Hz)')

hf = figure;
[stims,~,stimidx] = unique(labs,'stable');
cols = colorcube(10);
hold on;
for i = 1:numel(stims)
    idx = ismember(stimidx,i);
    x_mean = mean(data_naive(:,idx),1,'omitmissing');
    y_mean = mean(data_trained(:,idx),1,'omitmissing');
    x_std = std(data_naive(:,idx),0,1,'omitmissing')./sqrt(sum(idx));
    y_std = std(data_trained(:,idx),0,1,'omitmissing')./sqrt(sum(idx));
    
    % Plot errorbars
    errorbar(x_mean, y_mean, y_std, y_std, x_std, x_std, '.', ...
        'Color', cols(i,:), 'LineWidth', .5, 'HandleVisibility', 'off', ...
        'CapSize',1);
    
    % Plot scatter with size based on repnum
    scatter(x_mean, y_mean, 1 + 5./repnum(idx), ...
        'filled', 'CData', cols(i,:), 'MarkerFaceAlpha', 0.7, ...
        'DisplayName', stims{i});
end
maxrange = max([max(xlim),max(ylim)]);
minrange = min([min(xlim),min(ylim)]);
xlim([minrange maxrange]); ylim([minrange maxrange]);
plot([minrange maxrange],[minrange maxrange],'r-','DisplayName','x=y')
u = legendUnq; legend(u)
xline(0);yline(0);
hold off;
axis square
xlabel('naive iFR (Hz)');
ylabel('trained iFR (Hz)')
[~,p] = kstest2(data_naive(:),data_trained(:));
disp(['2-s KS Test of naive vs trained firing rates - p-val: ',num2str(p)])


v.dataFilter.stims_allowed = 'all stimuli';
v.dataFilter.interval = odor_interval;
v.dataFilter.subjectGroup = 'naïve';
out = v.plotUnitActivityMetricHead('method','avg intensity');
data_naive = cell2mat(out);
v.dataFilter.subjectGroup = 'trained';
out = v.plotUnitActivityMetricHead('method','avg intensity');
data_trained = cell2mat(out);

x = linspace(0,.6,1000);
figure; plotSplitViolin(gca, data_naive(:),data_trained(:), ...
    'show_median',true, ...
    'show_mean',true, ...
    'stats', true, ...
    'edges',x);
xlim([.5 1.5]); ylim([-.1 max(ylim)])
xticks([]); xlabel('naive | trained'); ylabel('iFR [Hz]')
set(gcf,'Color','w')
y = median(data_trained(:),'omitmissing') - median(data_naive(:),'omitmissing');
disp(['Trained - naive difference in median iFR: ' num2str(y), ' Hz'])
y = mean(data_trained(:),'omitmissing') - mean(data_naive(:),'omitmissing');
disp(['Trained - naive difference in mean iFR: ' num2str(y), ' Hz'])

% unit firing distributions (repetitions)
baseline_interval = [-22 -2];
odor_interval = dft.interval;
yrange = [-.3 .4];
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'trained';
v.dataFilter.trial_sorting = 'chronological';
v.dataFilter.stims_allowed = 'all stimuli';
[~,~, labs] = ModeSelector(v).extract; labs = labs{1};
% over reps
stims = unique(labs); nstims = numel(stims);
data_trained_allstim = [];
for i = 1:nstims
    v.dataFilter.stims_allowed = stims(i);
    v.dataFilter.interval = odor_interval;
    out = v.plotUnitActivityMetricHead('method','avg intensity');
    v.dataFilter.interval = baseline_interval;
    outbase = v.plotUnitActivityMetricHead('method','avg intensity');
    data_trained_allstim = [data_trained_allstim; cell2mat(out) - cell2mat(outbase)];
end
data = data_trained_allstim;
hf = figure;
subplot(121); b = prettyBoxplot(data,{'1','2','3','4','5'},'scatterSize',5,'plotLine',true);
ylabel('Cellwise delta iFR (Hz)')
ylim(yrange)
% Leu
v.dataFilter.stims_allowed = {'Leu'};
v.dataFilter.interval = odor_interval;
out = v.plotUnitActivityMetricHead('method','avg intensity');
v.dataFilter.interval = baseline_interval;
outbase = v.plotUnitActivityMetricHead('method','avg intensity');
data_trained_leu = cell2mat(out) - cell2mat(outbase);
subplot(122); b = prettyBoxplot(data_trained_leu,{'1','2','3','4','5'},'scatterSize',5,'plotLine',true);
ylim(yrange)
ylabel('Cellwise delta iFR (Hz)')
cfg.figSize = "medium";
cfg.aspRatioType = "square";
cfg.lineWidth = .5;
cfg.setFigure
cfg.saveFigure(gcf,['allgroups unit delta firing distribution repetitions raw'], saveType)
cfg = v.plotConfig;

% Friedman test: does firing change across repetitions? (trained)
fr_trained = statsUtils.friedman(data_trained_allstim);
fprintf('Friedman test (trained, all stimuli): chi2=%.2f, df=%d, p=%.4g\n', ...
    fr_trained.chi2, fr_trained.df, fr_trained.p);
fr_trained = statsUtils.friedman(data_trained_leu);
fprintf('Friedman test (trained, Leu): chi2=%.2f, df=%d, p=%.4g\n', ...
    fr_trained.chi2, fr_trained.df, fr_trained.p);

% repeat for naive group
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'naïve';
v.dataFilter.trial_sorting = 'chronological';
v.dataFilter.stims_allowed = 'all stimuli';
[~,~, labs_n] = ModeSelector(v).extract; labs_n = labs_n{1};
stims_n = unique(labs_n); nstims_n = numel(stims_n);
data_naive_allstim = [];
for i = 1:nstims_n
    v.dataFilter.stims_allowed = stims_n(i);
    v.dataFilter.interval = odor_interval;
    out = v.plotUnitActivityMetricHead('method','avg intensity');
    v.dataFilter.interval = baseline_interval;
    outbase = v.plotUnitActivityMetricHead('method','avg intensity');
    data_naive_allstim = [data_naive_allstim; cell2mat(out) - cell2mat(outbase)];
end
hf = figure;
subplot(121); b = prettyBoxplot(data_naive_allstim,{'1','2','3','4','5'},'scatterSize',5,'plotLine',true);
ylabel('Cellwise delta iFR (Hz)'); ylim(yrange); title('naive - all stimuli')
% Leu
v.dataFilter.stims_allowed = {'Leu'};
v.dataFilter.interval = odor_interval;
out = v.plotUnitActivityMetricHead('method','avg intensity');
v.dataFilter.interval = baseline_interval;
outbase = v.plotUnitActivityMetricHead('method','avg intensity');
data_naive_leu = cell2mat(out) - cell2mat(outbase);
subplot(122); b = prettyBoxplot(data_naive_leu,{'1','2','3','4','5'},'scatterSize',5,'plotLine',true);
ylim(yrange); ylabel('Cellwise delta iFR (Hz)'); title('naive - Leu')
cfg.figSize = "medium"; cfg.aspRatioType = "square"; cfg.lineWidth = .5;
cfg.setFigure
cfg.saveFigure(gcf,'naive unit delta firing distribution repetitions', saveType)
cfg = v.plotConfig;

fr_naive = statsUtils.friedman(data_naive_allstim);
fprintf('Friedman test (naive, all stimuli): chi2=%.2f, df=%d, p=%.4g\n', ...
    fr_naive.chi2, fr_naive.df, fr_naive.p);
fr_naive = statsUtils.friedman(data_naive_leu);
fprintf('Friedman test (naive, Leu): chi2=%.2f, df=%d, p=%.4g\n', ...
    fr_naive.chi2, fr_naive.df, fr_naive.p);
v.dataFilter = dft;

%% suppression score comparison across groups

yrange = [-4 4];

v.dataFilter = dft;
v.dataFilter.stims_allowed = {'Leu'};
v.dataFilter.subjectGroup = 'trained';
supp_trained = cell2mat(v.plotUnitActivityMetricHead('method','general suppression score'));
v.dataFilter.subjectGroup = 'naïve';
supp_naive = cell2mat(v.plotUnitActivityMetricHead('method','general suppression score'));
v.dataFilter = dft;

supp_stats = statsUtils.pairwiseMannWhitney({supp_naive, supp_trained}, {'naive','trained'}, true);
fprintf('Suppression score (naive vs trained): p_raw=%.4g, p_adj=%.4g, %s\n', ...
    supp_stats.p_raw, supp_stats.p_adj, statsUtils.pvalToStars(supp_stats.p_adj));
hf = figure; splitViolin({supp_naive, supp_trained},cfg,yrange);

cfg.figSize = 'small'; cfg.aspRatioType = 'square'; cfg.setFigure;
cfg.saveFigure(gcf,'suppression score Leu naive vs trained', saveType);
cfg = v.plotConfig;

v.dataFilter = dft;
v.dataFilter.stims_allowed = 'all stimuli';
v.dataFilter.subjectGroup = 'trained';
supp_trained = cell2mat(v.plotUnitActivityMetricHead('method','general suppression score'));
v.dataFilter.subjectGroup = 'naïve';
supp_naive = cell2mat(v.plotUnitActivityMetricHead('method','general suppression score'));
v.dataFilter = dft;

supp_stats = statsUtils.pairwiseMannWhitney({supp_naive, supp_trained}, {'naive','trained'}, true);
fprintf('Suppression score (naive vs trained): p_raw=%.4g, p_adj=%.4g, %s\n', ...
    supp_stats.p_raw, supp_stats.p_adj, statsUtils.pvalToStars(supp_stats.p_adj));
hf = figure; splitViolin({supp_naive, supp_trained},cfg,yrange);

cfg.figSize = 'small'; cfg.aspRatioType = 'square'; cfg.setFigure;
cfg.saveFigure(gcf,'suppression score allstims naive vs trained', saveType);
cfg = v.plotConfig;



function splitViolin(data,cfg,yrange)
    % split violin for naive (left) and trained (right)
    xcenter = 1;
    width = 0.25; % half-violin max width
    colors = cfg.c(1:2,:);
    
    hold on;
    for i = 1:2
        d = data{i};
        d = d(~isnan(d));
        if isempty(d), continue; end
        [f, xi] = ksdensity(d, 'NumPoints', 256);
        f = f / max(f) * width;
        if i == 1 % left (naive)
            X = [xcenter - f, xcenter*ones(1,numel(f))];
        else % right (trained)
            X = [xcenter + f, xcenter*ones(1,numel(f))];
        end
        Y = [xi, fliplr(xi)];
        patch(X, Y, colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.6);
        % % median line
        % med = median(d);
        % if i == 1
        %     plot([xcenter - width*0.9, xcenter], [med med], 'k', 'LineWidth', 1.5);
        % else
        %     plot([xcenter, xcenter + width*0.9], [med med], 'k', 'LineWidth', 1.5);
        % end
        % % jittered points
        % jitter = (rand(size(d)) - 0.5) * width * 0.6;
        % if i == 1
        %     xs = xcenter - 0.02 + jitter;
        % else
        %     xs = xcenter + 0.02 + jitter;
        % end
        % scatter(xs, d, 8, 'k', 'filled', 'MarkerFaceAlpha', 0.25);
    end
    xlim([xcenter-0.6, xcenter+0.6]);
    set(gca, 'XTick', xcenter, 'XTickLabel', {'naive    |    trained'});
    ylim(yrange); ylabel('Suppression (a.u.)');
    box off; axis square;
end



%
% % unit firing distributions
% baseline_interval = [-22 -2];
% odor_interval = dft.interval;
% yrange = [0 .2];
% v.dataFilter = dft;
% v.dataFilter.subjectGroup = 'naïve';
% v.dataFilter.trial_sorting = 'chronological';
% % over trials
% v.dataFilter.stims_allowed = 'all trials';
% hf = figure;
% out = v.plotTrialActivityMetricHead('method','avg intensity');
% ylim(yrange)
% % over repetitions
% v.dataFilter.stims_allowed = 'all stimuli';
% hf = figure;
% v.plotTrialActivityMetricHead('method','avg intensity', 'plotType', 'boxplot_repetitions')
% ylim(yrange)
% % over repetitions (Leu)
% v.dataFilter.stims_allowed = {'Leu'};
% hf = figure;
% v.plotTrialActivityMetricHead('method','avg intensity', 'plotType', 'boxplot_repetitions')
% ylim(yrange)
% 

%% general suppression score histogram
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'naïve';
c = v.plotUnitActivityMetricHead('method','general suppression score');
% histogram
figure; histogram(cell2mat(c),100, 'FaceColor','k','EdgeAlpha',0);
box off; axis square
xlabel('suppression [a.u]'); ylabel('histogram')
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
set(gcf, 'color', cfg.bgcol);

%% selectivity of tuning histogram
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

%% stability of tuning histogram
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


%% drift strength vs activity
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'naïve';
R = driftMetricsFigure(v,'plot',true);
drift = R.driftStrength3d;
close all
hf = figure;
nbins = 50;
xrange = [0 .4];
yrange = [-3 3];
xEdges = linspace(xrange(1),xrange(2),nbins+1);
yEdges = linspace(yrange(1),yrange(2),nbins+1);
% vs averaged activity across stimuli and reps
v.dataFilter.repetitions = 1:5;
out = v.plotUnitActivityMetricHead('metric','avg intensity');
firing = cell2mat(out); firing = mean(firing,[2,3],'omitmissing');
thisdrift = drift(:,:,1); thisdrift = mean(thisdrift,2,'omitmissing');
subplot(323); plotHeatmapAndIsoclines2(firing(:),thisdrift(:),xEdges,yEdges,1,0,1);
colormap(cfg.colormapName)
title('rep 1->2')
xlabel('Average iFR (Hz)'); ylabel('Delta iFR (Hz)')
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
[r,p] = corrcoef(firing(:),thisdrift(:));
disp(['Avg firing vs attenuation (rep 1>2): pears. r = ',num2str(r(2)),', p-val = ',num2str(p(2))])

thisdrift = drift(:,:,4); thisdrift = mean(thisdrift,2,'omitmissing');
subplot(324); plotHeatmapAndIsoclines2(firing(:),thisdrift(:),xEdges,yEdges,1,0,1);
colormap(cfg.colormapName)
title('rep 4->5')
xlabel('Average iFR (Hz)'); ylabel('Delta iFR (Hz)')
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
[r,p] = corrcoef(firing(:),thisdrift(:));
disp(['Avg firing vs attenuation (rep 4>5): pears. r = ',num2str(r(2)),', p-val = ',num2str(p(2))])
% vs stimulus-specific activity on the first trial of the pair
v.dataFilter.repetitions = 1;
out = v.plotUnitActivityMetricHead('metric','avg intensity');
xrange = [0 1];
yrange = [-15 15];
xEdges = linspace(xrange(1),xrange(2),nbins+1);
yEdges = linspace(yrange(1),yrange(2),nbins+1);
firing = cell2mat(out);
thisdrift = drift(:,:,1);
subplot(325); plotHeatmapAndIsoclines2(firing(:),thisdrift(:),xEdges,yEdges,1,0,1);
colormap(cfg.colormapName)
title('rep 1->2')
xlabel('iFR on rep 1 [Hz]'); ylabel('drift modulus')
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
v.dataFilter.repetitions = 4;
out = v.plotUnitActivityMetricHead('metric','avg intensity');
firing = cell2mat(out);
thisdrift = drift(:,:,v.dataFilter.repetitions);
subplot(326); plotHeatmapAndIsoclines2(firing(:),thisdrift(:),xEdges,yEdges,1,0,1);
colormap(cfg.colormapName)
title('rep 4->5')
xlabel('iFR on rep 4 [Hz]'); ylabel('drift modulus')
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);

cfg.figSize = "large";
cfg.aspRatioType = "square";
cfg.setFigure;
cfg.saveFigure(gcf,'trained drift vs activity', saveType)
cfg = v.plotConfig;