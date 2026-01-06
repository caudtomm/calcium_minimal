function [hf,data] = similarityDynamics(v,figs,windows,s)
% v : ExperimentViewer object
% figs : FigureSaver object
% windows : double [n,2] in seconds
% s : logical, saving option

dft = v.dataFilter;
cfg = v.plotConfig;


nwindows = height(windows);
out = cell(nwindows,1);
for i = 1:nwindows
    v.dataFilter.interval = windows(i,:);
    [hf,out{i}, subject_groups, stim_groups] = plotRepetitionDistances(v,'correlation'); % outputs 2 figures
    figs.title = ['Repetitions: sec', num2str(windows(i,1)), '-', num2str(windows(i,2))];
    % if s; figs.append(hf); end
    close(hf)
end
v.dataFilter = dft;

% extract relevant data
nplots = numel(out{1});
data = nan(nwindows,nplots,3); % dim3: [mean, sem, std]
for i_w = 1:nwindows
    for i_p = 1:nplots
        thisplot = out{i_w}{i_p};
        if isempty(thisplot); continue; end
        data(i_w,i_p,1) = mean(thisplot.data(:),'omitmissing'); % mean
        data(i_w,i_p,2) = std(thisplot.data(:),[],'omitmissing')/numel(thisplot.data(:)); % sem
        data(i_w,i_p,3) = std(thisplot.data(:),[],'omitmissing'); % std
    end
end

%% plotting

% make sure stim_groups entries are usable as labels
for i = 1:numel(stim_groups)
    if iscell(stim_groups{i})
        stim_groups{i} = strjoin(stim_groups{i}, '');
    end
end

% a couple useful vars
n_stimgroups = numel(stim_groups);
n_subject_groups = numel(subject_groups);
idx_by_stimgroup = repmat(1:n_stimgroups,1,n_subject_groups);
idx_by_subjectgroup = repelem(1:n_subject_groups,1,n_stimgroups);

% all together
hf = figure;
b = plotCurves(data,windows,cfg);

% one subplot per subject-group
hf = figure;
for i = 1:n_subject_groups
    subplot(n_subject_groups,1,i)
    b = plotCurves(data(:, idx_by_subjectgroup==i, :),windows,cfg);
    if i==1; legend(b,stim_groups); end
    title(subject_groups{i})
end

% one subplot per stim-group
hf = figure;
for i = 1:n_stimgroups
    subplot(n_stimgroups,1,i)
    b = plotCurves(data(:, idx_by_stimgroup==i, :),windows,cfg);
    if i==1; legend(b,subject_groups); end
    title(stim_groups{i})
end


%% plot all familiar in trained vs all novel in trained vs all novel in naive

t = mean(windows,2);
[~,tstart] = min(abs(t)); % t0 = stim onset
% tstart = tstart + 1; % <- actually, the immediately following frame

idx = find(idx_by_stimgroup==6 & idx_by_subjectgroup==2); % all familiar & trained
p1 = getFitParams(out, idx, t, tstart);

idx = find(idx_by_stimgroup==7 & idx_by_subjectgroup==2); % all novel & trained
p2 = getFitParams(out, idx, t, tstart);

idx = find(idx_by_stimgroup==7 & idx_by_subjectgroup==1); % all novel & naive
p3 = getFitParams(out, idx, t, tstart);

hf = figure;
set(gcf, 'color', cfg.bgcol); 

subplot(161); hold on; clear b % # ----------------------

% trained and familiar
b(1) = scatter(t, p1.y_avg,'b','filled');
errorbar(t, p1.y_avg,p1.err, ...
    'vertical', 'LineStyle', 'none','Color','b');
plot(p1.avgfit.t, p1.avgfit.vals,'b--')

% trained and novel
b(2) = scatter(t, p2.y_avg,'r','filled');
errorbar(t, p2.y_avg,p2.err, ...
    'vertical', 'LineStyle', 'none','Color','r');
plot(p2.avgfit.t, p2.avgfit.vals,'r--')

b(3) = scatter(t, p3.y_avg,'g','filled');
errorbar(t, p3.y_avg,p3.err, ...
    'vertical', 'LineStyle', 'none','Color','g');
plot(p3.avgfit.t, p3.avgfit.vals,'g--')

% naive and novel
axis tight
legend(b,{'familiar/trained','novel/trained','novel/naive'})
xlabel('Time from stim. onset [s]');
ylabel('Avg. intertrial similarity (same odor) + SEM');
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
% # ----------------------
data = [];
for i = 1:3
    thisp = sprintf('p%s',num2str(i));
    p.amp{i} = cellfun(@(x) x.p(1), eval([thisp,'.win_fits']));
    p.tau{i} = cellfun(@(x) 1./x.p(2), eval([thisp,'.win_fits']));
    p.offset{i} = cellfun(@(x) x.p(3), eval([thisp,'.win_fits']));
    p.half_t{i} = cellfun(@(x) x.half_t, eval([thisp,'.win_fits']));
    p.t2max{i} = cellfun(@(x) x.t2max, eval([thisp,'.win_fits']));

    data = [data ; [p.amp{i},p.tau{i},p.offset{i},p.half_t{i},p.t2max{i}]];
end
g = [repmat({'trained-familiar'},numel(p.amp{1}),1) ; ...
    repmat({'trained-novel'},numel(p.amp{2}),1) ; ...
    repmat({'naive'},numel(p.amp{3}),1)];
labs = {'amplitude','tau','offset','half t','T'};
grouplabs = {'trained-familiar','trained-novel','naive'};

avg_dt = [[p1.avgfit.p,p1.avgfit.half_t];
          [p2.avgfit.p,p2.avgfit.half_t];
          [p3.avgfit.p,p3.avgfit.half_t]];

for i = 2:5
    subplot(1,6,i)
    boxplot(data(:,i-1),g);
    hold on
    scatter(1:3,avg_dt(:,i-1),50,'red','filled')
    box off
    ylabel(labs{i-1});
    set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
end
% 
% subplot(162)
% violin(p.amp)
% box off
% xticklabels(grouplabs)
% ylabel(labs{1});
% set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
% 
% subplot(163)
% violin(p.tau)
% box off
% xticklabels(grouplabs)
% ylabel(labs{2});
% set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
% 
% subplot(164)
% violin(p.offset)
% box off
% xticklabels(grouplabs)
% ylabel(labs{3});
% set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
% 
% subplot(165)
% violin(p.half_t)
% box off
% xticklabels(grouplabs)
% ylabel(labs{4});
% set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);

subplot(166)
bar([p.t2max{1}(1),p.t2max{2}(1),p.t2max{3}(1)])
box off
xticklabels(grouplabs)
ylabel(labs{5});
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);




end

function b = plotCurves(data,windows,cfg)
% plot onto provided axes

% Get dimensions
[n_time, n_vars, ~] = size(data);
t = mean(windows,2)'; % time axis


% Plot with SEM whiskers
hold on;
b = gobjects(n_vars,1);
for i = 1:n_vars
    mu = data(:, i, 1);
    sem = data(:, i, 2);
    
    % Shaded error bars (optional, looks cleaner)
    fill([t fliplr(t)], [mu - sem; flipud(mu + sem)]', ...
         cfg.c(i,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % grey shade
    
    % Plot mean line
    b(i) = plot(t, mu, 'LineWidth', cfg.lineWidth, 'Color', cfg.c(i,:));
end
axis tight
xlabel('Time from stim. onset');
ylabel('Avg. intertrial similarity (same odor)');
title('Mean ± SEM');
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
set(gcf, 'color', cfg.bgcol); 
end


function params = getFitParams(out, idx, t, tstart)

% -----------------
nwindows = numel(out);
[ncomparisons,ndatapoints] = size(out{1}{idx}.data);
G = nan(ncomparisons,ndatapoints,nwindows); % to permute later
for i = 1:nwindows; G(:,:,i) = out{i}{idx}.data; end
G = permute(G,[3,2,1]); % [windows, datapoints(odors*subjects), comparisons(intertrial corrs)]
% -----------------
g = mean(G,3,"omitmissing");

y_avg = mean(g,2,'omitmissing');
err = std(g,[],2,'omitmissing')./ndatapoints; % sem
[~,tend] = max(y_avg);
tend = tend+1;
data_t = t(tstart:tend);

win_fits = cell(ndatapoints,1);
for i = 1:ndatapoints
    [win_fits{i},model] = fitModel(data_t,g(:,i),t,tstart,tend);
end

params.y_avg = y_avg;
params.err = err;
params.model = model;
params.t = data_t;
params.avgfit = fitModel(data_t,y_avg,t,tstart,tend);
params.win_fits = win_fits;

end


function [fit,model] = fitModel(data_t,y,t,tstart,tend)
[p, ~, model] = fitExpSaturation(data_t, y(tstart:tend), 0);
t_fit = t(tstart):.01:t(tend);
yfit = model(p,t_fit);
[~,half_t] = min(abs(yfit-(min(yfit)+(max(yfit)-min(yfit))/2))); half_t = t_fit(half_t);
t2max = t(tend);

fit.p = p;
fit.t = t_fit;
fit.vals = yfit;
fit.half_t = half_t;
fit.t2max = t2max;

end