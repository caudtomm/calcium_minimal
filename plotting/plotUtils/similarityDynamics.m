function [hf,data] = similarityDynamics(v,figs,windows,s)
% v : ExperimentViewer object
% figs : FigureSaver object
% windows : double [n,2] in seconds
% s : logical, saving option

cfg = v.plotConfig;


nwindows = height(windows);
out = cell(nwindows,1);
for i = 1:nwindows
    [hf,out{i}, subject_groups, stim_groups] = v.plotRepetitionDistances(windows(i,:),'correlation'); % outputs 2 figures
    figs.title = ['Repetitions: sec', num2str(windows(i,1)), '-', num2str(windows(i,2))];
    if s; figs.append(hf); end
    close(hf)
end

% extract relevant data
nplots = numel(out{1});
data = nan(nwindows,nplots,3); % dim3: [mean, sem, std]
for i_w = 1:nwindows
    for i_p = 1:nplots
        thisplot = out{i_w}{i_p};
        if isempty(thisplot); continue; end
        data(i_w,i_p,1) = mean(thisplot.data,'omitmissing'); % mean
        data(i_w,i_p,2) = std(thisplot.data,[],'omitmissing')/numel(thisplot.data); % sem
        data(i_w,i_p,3) = std(thisplot.data,[],'omitmissing'); % std
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
tstart = tstart + 1; % <- actually, the immediately following frame

idx = find(idx_by_stimgroup==6 & idx_by_subjectgroup==2); % all familiar & trained
y1 = data(:,idx,1);
err_y1 = data(:,idx,3);
[~,tend1] = max(y1);
% tend1 = tend1+1;
[params1, yfit1, ~] = fitExpSaturation(t(tstart:tend1), y1(tstart:tend1), 0);

idx = find(idx_by_stimgroup==7 & idx_by_subjectgroup==2); % all novel & trained
y2 = data(:,idx,1);
err_y2 = data(:,idx,3);
[~,tend2] = max(y2);
% tend2 = tend2+1;
[params2, yfit2, ~] = fitExpSaturation(t(tstart:tend2), y2(tstart:tend2), 0);

idx = find(idx_by_stimgroup==7 & idx_by_subjectgroup==1); % all novel & naive
y3 = data(:,idx,1);
err_y3 = data(:,idx,3);
[~,tend3] = max(y3);
% tend3 = tend3+1;
[params3, yfit3, ~] = fitExpSaturation(t(tstart:tend3), y3(tstart:tend3), 0);

hf = figure;
set(gcf, 'color', cfg.bgcol); 

subplot(121); hold on; clear b
b(1) = scatter(t(tstart:tend1), y1(tstart:tend1),'b','filled');
errorbar(t(tstart:tend1), y1(tstart:tend1),err_y1(tstart:tend1), ...
    'vertical', 'LineStyle', 'none','Color','b');
b(2) = scatter(t(tstart:tend2), y2(tstart:tend2),'r','filled');
errorbar(t(tstart:tend2), y2(tstart:tend2),err_y2(tstart:tend2), ...
    'vertical', 'LineStyle', 'none','Color','r');
b(3) = scatter(t(tstart:tend3), y3(tstart:tend3),'g','filled');
errorbar(t(tstart:tend3), y3(tstart:tend3),err_y3(tstart:tend3), ...
    'vertical', 'LineStyle', 'none','Color','g');
plot(t(tstart:tend1), yfit1,'b--')
plot(t(tstart:tend2), yfit2,'r--')
plot(t(tstart:tend3), yfit3,'g--')
axis tight
legend(b,{'familiar/trained','novel/trained','novel/naive'})
xlabel('Time from stim. onset [s]');
ylabel('Avg. intertrial similarity (same odor) + STD');
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);

subplot(122)
bar([params1;params2;params3]')
xticklabels({'amplitude','1/tau','offset'})
box off
ylabel('Fitted value');
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