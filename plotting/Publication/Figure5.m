dbstop if error

s = true; % save figures to files?
savepath = 'bin5';
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
 


%% FIGURE 5

%% dPCA decomposition
v.plotConfig.custom.crange = [.2 1];

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
v.dataFilter.mode_OI = 'stimulus';
hf = figure;
outMat_up = v.plotDistancesHead;
v.dataFilter.mode_OI = 'novelty';
v.dataFilter.interval = [0 3];
hf = figure;
c = v.plotTrialActivityMetricHead('method','avg intensity', 'plotType', 'boxplot_repetitions');
v.dataFilter.interval = [-1 3];
[~,events,labs] = ModeSelector(v).extract;
events = cellfun(@squeeze,events,'UniformOutput',false);
events = cell2mat(permute(events,[2,3,1]));
t = linspace(v.dataFilter.interval(1),v.dataFilter.interval(2),height(events));
labs = labs{1};
hf = figure;
imagesc(t,1:numel(labs),mean(events,3,'omitmissing')')
yticks(1:numel(labs)); yticklabels(labs);
colormap(cfg.colormapName)
%
v.dataFilter = dft;


%% stimulus dPC weights
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


%% stimulus dPC weights
v.dataFilter = dft;
v.dataFilter.mode_name = 'dpca';
v.dataFilter.mode_OI = 'stimulus';
v.dataFilter.mode_method = 'mode_values';
v.dataFilter.subjectGroup = 'all';
v.dataFilter.mode_file = 'dpca_shuffle_all.mat';
m = ModeSelector(v).extract;
c = [];
for i= 1:numel(m.coeffs)
    thisdata = m.fullout{i}.V(:,m.parseModeOI(i));
    thisdata = sum(thisdata.^2);
    c = [c; thisdata(:)];
end
%histogram
figure; histogram(c, 100, 'FaceColor','k','EdgeAlpha',0);
box off; axis square
axis tight
xlabel('Sum of sq. weights'); ylabel('Histogram')
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
set(gcf, 'color', cfg.bgcol);
% 
c = []; % [w1 N; w2 N; ...]
for i= 1:numel(m.coeffs)
    thisdata = m.fullout{i}.V(:,m.parseModeOI(i));
    c = [c; thisdata(:), ones(numel(thisdata),1)*v.filtered_traces{i}.N];
end
%histogram
figure; histogram(abs(c(:,1)).*c(:,2), 100, 'FaceColor','k','EdgeAlpha',0);
box off; axis square
yscale log
xscale log
hold on; plot([1 1],[min(ylim) max(ylim)],"Color",'k') % null hypothesis: all |weights| equal
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
box off; axis square tight
xlim([-1 1])
xlabel('w similarity'); ylabel('histogram')
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
set(gcf, 'color', cfg.bgcol);


v.dataFilter = dft;
v.dataFilter.mode_name = 'dpca';
v.dataFilter.mode_OI = 'stimulus';
v.dataFilter.mode_method = 'mode_values';
v.dataFilter.subjectGroup = 'all';
v.dataFilter.mode_file = 'dpca_shuffle_all.mat';
m = ModeSelector(v).extract;
v.dataFilter.mode_OI = 'novelty';
n = ModeSelector(v).extract;
c = [];
for i= 1:nsubjects
    thisstimw = m.coeffs{i}(:,m.parseModeOI(i));
    thisnovw = n.coeffs{i}(:,n.parseModeOI(i));
    for j = 1:width(thisstimw)
        thisdata = [1-pdist([thisstimw(:,j), thisnovw]','cosine')];
        c = [c; thisdata];
    end
end
figure;
subplot(211); histogram(c,100,'FaceColor','k','EdgeAlpha',0);
box off; axis square tight
xlim([-1 1])
xlabel('w similarity'); ylabel('histogram')
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
set(gcf, 'color', cfg.bgcol);
v.dataFilter = dft;
v.dataFilter.mode_name = 'dpca';
v.dataFilter.mode_OI = 'novelty';
v.dataFilter.mode_method = 'mode_values';
v.dataFilter.subjectGroup = 'all';
v.dataFilter.mode_file = 'dpca_shuffle_all.mat';
m = ModeSelector(v).extract;
c = []; % [w1 N; w2 N; ...]
for i= 1:numel(m.coeffs)
    thisdata = m.coeffs{i}(:,m.parseModeOI(i));
    c = [c; thisdata(:), ones(numel(thisdata),1)*v.filtered_traces{i}.N];
end
% histogram
y = abs(c(:,1)).*c(:,2);
subplot(212); histogram(y, 100, 'FaceColor','k','EdgeAlpha',0);
[~,p,~] = kstest(y-1);
disp(['1-sample KS test - p-val: ',num2str(p)])
box off; axis square
yscale log
xscale log
hold on; plot([1 1],[min(ylim) max(ylim)],"Color",'k') % null hypothesis: all |weights| equal
plot(mean(y)*[1 1],[min(ylim) max(ylim)],"Color",'r') % average norm w
axis tight
xlabel('|w|*N'); ylabel('histogram')
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
set(gcf, 'color', cfg.bgcol);
cfg.figSize = 'small';
cfg.aspRatioType = 'tall';
cfg.setFigure;
cfg.saveFigure(gcf,'naive novelty weights', saveType)



v.dataFilter = dft;
v.dataFilter.mode_name = 'dpca';
v.dataFilter.mode_OI = 'stimulus';
v.dataFilter.mode_method = 'mode_values';
v.dataFilter.subjectGroup = 'all';
v.dataFilter.mode_file = 'dpca_shuffle_all.mat';
m = ModeSelector(v).extract;
c = []; % []
for i= 1:numel(m.coeffs)
    thisdata = m.coeffs{i}(:,m.parseModeOI(i));
    thisdata = thisdata.^2;
    sumw = sum(thisdata(:));
    stdw = std(thisdata(:));
    c = [c; (stdw^2)];
end
hf = figure;
y = c(:); histogram(y, 100, 'FaceColor','k','EdgeAlpha',0);
% [~,p,~] = kstest(y-1);
% disp(['1-sample KS test - p-val: ',num2str(p)])
box off; axis square tight
xlabel('Var(w)'); ylabel('histogram')
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
set(gcf, 'color', cfg.bgcol);
cfg.figSize = 'tiny';
cfg.aspRatioType = 'square';
cfg.setFigure;
cfg.saveFigure(gcf,'naive novelty weights', saveType)


%% variance explained by each PC/dPC
v.dataFilter = dft;
v.dataFilter.mode_name = 'dpca';
v.dataFilter.mode_OI = 'all_stimulus';
v.dataFilter.mode_method = 'mode_values';
v.dataFilter.subjectGroup = 'naïve';
v.dataFilter.mode_file = 'dpca_naive.mat';
m = ModeSelector(v).extract;
nsubjects = numel(m.fullout);
data = nan(3,20,nsubjects);
n = zeros(nsubjects,1);
for i=1:nsubjects
    thisdata = m.fullout{i}.explVar;
    data(1,:,i) = thisdata.cumulativePCA;
    
    idx = m.parseModeOI(i);
    idx_novelty = find(idx,1);
    idx(idx_novelty)=false;
    n(i) = sum(idx);

    data(2,1:sum(idx),i) = cumsum(thisdata.componentVar(idx));
    data(3,1,i) = thisdata.componentVar(idx_novelty);
    
end
disp(['# of identity dPCs: ',num2str(min(n)),'-',num2str(max(n))])
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


%% novelty mode weights
v.dataFilter = dft;
v.dataFilter.mode_name = 'dpca';
v.dataFilter.mode_OI = 'novelty';
v.dataFilter.mode_method = 'mode_values';
v.dataFilter.subjectGroup = 'trained';
v.dataFilter.mode_file = 'dpca_trained.mat';
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
%% novelty w correlation with general suppression
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'trained';
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
%% novelty w correlation with odor-specific suppression
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
%% novelty w correlation with population intensity
v.dataFilter = dft;
v.dataFilter.subjectGroup = 'trained';
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
