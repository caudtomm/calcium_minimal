function [hf, data] = similarityDynamics(v, figs, windows, s, subject_groups, stim_groups, err_type, fit_stim_keys, fit_subject_keys)
% v                : ExperimentViewer
% figs             : FigureSaver
% windows          : [n×2] double, time windows in seconds
% s                : logical, saving flag (reserved)
% subject_groups   : cell of subject group names
% stim_groups      : cell of stim group specifiers (strings or cell-arrays of odor names)
% err_type         : 'sem' or 'std' — error ribbon type for time-course plots
% fit_stim_keys    : cell of stim group specifiers to fit ([] = all stim groups)
%                    each element is a string or cell-array of odor names, matching an entry in stim_groups
% fit_subject_keys : cell of subject group names for fitting ({}  = all subject groups)
%                    one figure is generated per stim group; all subject groups appear in each figure
arguments
    v ExperimentViewer
    figs
    windows double
    s logical = false
    subject_groups cell = {'naïve', 'trained'}
    stim_groups cell = {{'Arg','Ala','His'}, {'Trp','Ser','Leu'}, 'all stimuli'}
    err_type string {mustBeMember(err_type, {'sem','std'})} = 'sem'
    fit_stim_keys = []
    fit_subject_keys cell = {}
end

dft = v.dataFilter;
cfg = v.plotConfig;

%% Collect data across time windows

nwindows = height(windows);
out = cell(nwindows, 1);
for i = 1:nwindows
    v.dataFilter.interval = windows(i,:);
    [hf, out{i}, subject_groups, stim_groups] = plotRepetitionDistances(v, 'correlation', subject_groups, stim_groups);
    % if s; figs.append(hf); end
    close(hf)
end
v.dataFilter = dft;

%% Index maps and display labels

stim_labels      = cellfun(@stimGroupLabel, stim_groups, 'UniformOutput', false);
n_stim_groups    = numel(stim_groups);
n_subject_groups = numel(subject_groups);
nplots           = n_stim_groups * n_subject_groups;

% plot ordering matches plotRepetitionDistances: outer loop = subject, inner = stim
idx_by_subjectgroup = repelem(1:n_subject_groups, 1, n_stim_groups);
idx_by_stimgroup    = repmat(1:n_stim_groups,    1, n_subject_groups);

%% Extract summary statistics per (window, plot)

err_slice = 2 + strcmp(err_type, 'std'); % 2 = sem, 3 = std
data = nan(nwindows, nplots, 3); % dim3: [mean, sem, std]
for i_w = 1:nwindows
    for i_p = 1:nplots
        thisplot = out{i_w}{i_p};
        if isempty(thisplot); continue; end
        vals  = thisplot.data(:);
        n_obs = sum(~isnan(vals));
        data(i_w, i_p, 1) = mean(vals, 'omitmissing');
        data(i_w, i_p, 2) = std(vals, [], 'omitmissing') / sqrt(n_obs);
        data(i_w, i_p, 3) = std(vals, [], 'omitmissing');
    end
end

%% Time-course curves: all conditions together

hf_all = figure;
plotCurves(data, windows, cfg, err_slice);

%% Time-course curves: one subplot per subject group

hf_bysubj = figure;
for i = 1:n_subject_groups
    subplot(n_subject_groups, 1, i)
    b = plotCurves(data(:, idx_by_subjectgroup==i, :), windows, cfg, err_slice);
    if i == 1; legend(b, stim_labels); end
    title(subject_groups{i})
end

%% Time-course curves: one subplot per stim group

hf_bystim = figure;
for i = 1:n_stim_groups
    subplot(n_stim_groups, 1, i)
    b = plotCurves(data(:, idx_by_stimgroup==i, :), windows, cfg, err_slice);
    if i == 1; legend(b, subject_groups); end
    title(stim_labels{i})
end

%% Per-time-bin pairwise group comparisons (FDR-adjusted Wilcoxon rank-sum)

t_centers = mean(windows, 2);
res = cell(n_stim_groups, 1);
for i_sg = 1:n_stim_groups
    pair_labels = {};
    for a = 1:n_subject_groups
        for jg = a+1:n_subject_groups
            pair_labels{end+1} = [subject_groups{a}, ' vs ', subject_groups{jg}]; %#ok
        end
    end
    npairs   = numel(pair_labels);
    p_matrix = nan(nwindows, npairs);

    for i_w = 1:nwindows
        grp_data = cell(n_subject_groups, 1);
        for g = 1:n_subject_groups
            i_p = find(idx_by_stimgroup==i_sg & idx_by_subjectgroup==g);
            if ~isempty(i_p) && ~isempty(out{i_w}{i_p(1)})
                grp_data{g} = out{i_w}{i_p(1)}.data(:);
            end
        end
        pidx = 0;
        for a = 1:n_subject_groups
            for jg = a+1:n_subject_groups
                pidx = pidx + 1;
                if ~isempty(grp_data{a}) && ~isempty(grp_data{jg})
                    p_matrix(i_w, pidx) = ranksum(grp_data{a}, grp_data{jg});
                end
            end
        end
    end

    p_adj_matrix = p_matrix;
    for pidx = 1:npairs
        p_adj_matrix(:, pidx) = statsUtils.fdr(p_matrix(:, pidx));
    end

    T = table(t_centers, 'VariableNames', {'t_center'});
    for pidx = 1:npairs
        T.(matlab.lang.makeValidName(pair_labels{pidx})) = p_adj_matrix(:, pidx);
    end
    fprintf('\n=== Per-bin group comparison (stim group: %s) — FDR-adjusted p-values ===\n', stim_labels{i_sg});
    disp(T);
    res{i_sg} = T;
end

%% P-value over time: one figure per subject-group pair

n_pairs    = n_subject_groups * (n_subject_groups - 1) / 2;
hf_pvalues = gobjects(n_pairs, 1);
ip         = 0;
if n_subject_groups >= 2
    for a = 1:n_subject_groups
        for jg = a+1:n_subject_groups
            pair_col   = matlab.lang.makeValidName([subject_groups{a}, ' vs ', subject_groups{jg}]);
            pair_title = [subject_groups{a}, ' vs ', subject_groups{jg}];
            ip = ip + 1;
            hf_pvalues(ip) = figure; hold on;
            b_pv = gobjects(n_stim_groups, 1);
            for i_sg = 1:n_stim_groups
                if ismember(pair_col, res{i_sg}.Properties.VariableNames)
                    b_pv(i_sg) = plot(res{i_sg}.t_center, res{i_sg}.(pair_col), 'Color', cfg.c(i_sg,:));
                end
            end
            plot(xlim, [.05 .05], 'r--')
            yscale log
            xlabel('Time from stim. onset (s)')
            ylabel('P-value (FDR-adjusted)')
            title(pair_title)
            legend(b_pv, stim_labels)
            cfg.figSize = 'small'; cfg.aspRatioType = 'wide'; cfg.setFigure;
            cfg.saveFigure(gcf, ['similarity dynamics - pvalues - ', pair_title], 'vector')
        end
    end
end

%% Exponential fitting: one figure per stim group, comparing all specified subject groups
% [] / {} expand to all groups; wrap cell-array stim specifiers in an outer cell,
% e.g. fit_stim_keys = {{'Trp','Ser','Leu'}, 'all stimuli'}

if isempty(fit_stim_keys)
    fit_stim_list = stim_groups;
else
    fit_stim_list = fit_stim_keys;
end
if isempty(fit_subject_keys)
    fit_subj_list = subject_groups;
else
    fit_subj_list = fit_subject_keys;
end

n_fit_stims  = numel(fit_stim_list);
n_fit_groups = numel(fit_subj_list);
hf_fits      = gobjects(n_fit_stims, 1);
mode = 'exponential';
t    = mean(windows, 2);
[~, tstart] = min(abs(t));

for i_fs = 1:n_fit_stims
    i_stim         = findGroupIdx(stim_groups, fit_stim_list{i_fs});
    stim_label_fit = stimGroupLabel(fit_stim_list{i_fs});

    fit_params = cell(n_fit_groups, 1);
    for k = 1:n_fit_groups
        i_subj = find(strcmp(subject_groups, fit_subj_list{k}));
        if isempty(i_subj)
            warning('Subject group not found for fitting: %s', fit_subj_list{k});
            continue
        end
        i_plot = find(idx_by_stimgroup==i_stim & idx_by_subjectgroup==i_subj(1));
        if ~isempty(i_plot)
            fit_params{k} = getFitParams(out, i_plot(1), t, tstart, mode);
        end
    end

    hf_fits(i_fs) = figure;
    set(gcf, 'color', cfg.bgcol);
    sgtitle(stim_label_fit, 'Color', cfg.textcol)

    subplot(1,6,1); hold on;
    b_fit = gobjects(n_fit_groups, 1);
    for k = 1:n_fit_groups
        if isempty(fit_params{k}); continue; end
        pk = fit_params{k};
        b_fit(k) = scatter(t, pk.y_avg, 'filled', ...
            'MarkerFaceColor', cfg.c(k,:), 'MarkerEdgeColor', cfg.c(k,:));
        errorbar(t, pk.y_avg, pk.err, 'vertical', 'LineStyle', 'none', 'Color', cfg.c(k,:));
        plot(pk.avgfit.t, pk.avgfit.vals, '--', 'Color', cfg.c(k,:))
    end
    axis tight; xlim([-1 3])
    legend(b_fit, fit_subj_list)
    xlabel('Time from stim. onset [s]')
    ylabel('Avg. intertrial similarity ± SEM')
    set(gca, 'color', cfg.bgcol, 'XColor', cfg.axcol, 'YColor', cfg.axcol, 'ZColor', cfg.axcol);

    params_data      = [];
    group_labels_fit = {};
    avg_dt           = nan(n_fit_groups, 4); % [amp, tau, offset, half_t] per group
    for k = 1:n_fit_groups
        if isempty(fit_params{k}); continue; end
        pk   = fit_params{k};
        amp  = cellfun(@(x) x.p(1),   pk.win_fits);
        tau  = cellfun(@(x) 1/x.p(2), pk.win_fits);
        offs = cellfun(@(x) x.p(3),   pk.win_fits);
        ht   = cellfun(@(x) x.half_t, pk.win_fits);
        params_data      = [params_data;      amp(:), tau(:), offs(:), ht(:)]; %#ok
        group_labels_fit = [group_labels_fit; repmat(fit_subj_list(k), numel(amp), 1)]; %#ok
        avg_dt(k,:)      = [pk.avgfit.p(1), 1/pk.avgfit.p(2), pk.avgfit.p(3), pk.avgfit.half_t];
    end

    param_labels = {'amplitude','tau','offset','half t'};
    for i = 1:4
        subplot(1, 6, i+1)
        boxplot(params_data(:,i), group_labels_fit);
        hold on
        scatter(1:n_fit_groups, avg_dt(:,i), 50, 'red', 'filled')
        box off
        ylabel(param_labels{i})
        set(gca, 'color', cfg.bgcol, 'XColor', cfg.axcol, 'YColor', cfg.axcol, 'ZColor', cfg.axcol);
    end

    subplot(1,6,6)
    t2max_vals = nan(n_fit_groups, 1);
    for k = 1:n_fit_groups
        if ~isempty(fit_params{k}); t2max_vals(k) = fit_params{k}.avgfit.t2max; end
    end
    bar(t2max_vals)
    box off
    xticklabels(fit_subj_list)
    ylabel('Time to max [s]')
    set(gca, 'color', cfg.bgcol, 'XColor', cfg.axcol, 'YColor', cfg.axcol, 'ZColor', cfg.axcol);
end

hf = struct( ...
    'timecourse_all',    hf_all, ...
    'timecourse_bysubj', hf_bysubj, ...
    'timecourse_bystim', hf_bystim, ...
    'pvalues',           hf_pvalues, ...
    'fits',              hf_fits);

end


function b = plotCurves(data, windows, cfg, err_slice)
[~, n_vars, ~] = size(data);
t = mean(windows, 2)';

hold on;
b = gobjects(n_vars, 1);
for i = 1:n_vars
    mu  = data(:, i, 1);
    err = data(:, i, err_slice);
    fill([t, fliplr(t)], [mu-err; flipud(mu+err)]', cfg.c(i,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    b(i) = plot(t, mu, 'LineWidth', cfg.lineWidth, 'Color', cfg.c(i,:));
end
axis tight
xlabel('Time from stim. onset')
ylabel('Avg. intertrial similarity (same odor)')
set(gca, 'color', cfg.bgcol, 'XColor', cfg.axcol, 'YColor', cfg.axcol, 'ZColor', cfg.axcol);
set(gcf, 'color', cfg.bgcol);
end


function params = getFitParams(out, idx, t, tstart, mode)
nwindows = numel(out);
[ncomparisons, ndatapoints] = size(out{1}{idx}.data);
G = nan(ncomparisons, ndatapoints, nwindows);
for i = 1:nwindows; G(:,:,i) = out{i}{idx}.data; end
G = permute(G, [3,2,1]); % [windows, datapoints, comparisons]

g     = mean(G, 3, 'omitmissing');
y_avg = mean(g, 2, 'omitmissing');
err   = std(g, [], 2, 'omitmissing') / sqrt(ndatapoints); % sem across subjects/odors

[~, tend] = max(y_avg);
tend   = tend + 1;
data_t = t(tstart:tend);

win_fits = cell(ndatapoints, 1);
for i = 1:ndatapoints
    [win_fits{i}, model] = fitModel(data_t, g(:,i), t, tstart, tend, mode);
end

params.y_avg    = y_avg;
params.err      = err;
params.model    = model;
params.t        = data_t;
params.avgfit   = fitModel(data_t, y_avg, t, tstart, tend, mode);
params.win_fits = win_fits;
end


function [fit, model] = fitModel(data_t, y, t, tstart, tend, mode)
switch mode
    case 'exponential'
        [p, ~, model] = fitExpSaturation(data_t, y(tstart:tend), 0);
    case '2-point line'
        [p, ~, model] = lineThrough2Points(data_t, y(tstart:tend));
    otherwise
        error('Model class not recognized.')
end

t_fit = t(tstart):.01:t(tend);
yfit  = model(p, t_fit);
[~, half_t_idx] = min(abs(yfit - (min(yfit) + (max(yfit)-min(yfit))/2)));

fit.p      = p;
fit.t      = t_fit;
fit.vals   = yfit;
fit.half_t = t_fit(half_t_idx);
fit.t2max  = t(tend);
end


function [params, y_fit, model] = lineThrough2Points(t, y)
model  = @(p, x) p(1)*x + p(2);
p1     = (y(end) - y(1)) / (t(end) - t(1));
p2     = y(1) - p1*t(1);
params = [p1, p2];
y_fit  = model(params, t);
end


function label = stimGroupLabel(group)
    if iscell(group)
        label = strjoin(group, ', ');
    else
        label = group;
    end
end


function idx = findGroupIdx(groups, key)
    idx = find(cellfun(@(g) isequal(g, key), groups));
    if isempty(idx); error('Stim group key not found.'); end
    idx = idx(1);
end
