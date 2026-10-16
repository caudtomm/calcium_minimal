function h = pairedScatter(ax, data, opts)
% pairedScatter  Jittered scatter with paired connecting lines per column.
%
%   h = pairedScatter(ax, data)
%   h = pairedScatter(ax, data, Name, Value, ...)
%
%   Each row of data is one variable (subject, cell, …); each column is one
%   measurement/condition. Points in the same row are connected by a line.
%
%   Inputs:
%     ax   - axes handle, or [] for gca
%     data - [variables x measurements] double matrix
%
%   Name-Value Options:
%     xpos        - [1 x M] x-axis positions              (default: 1,2,…,M)
%     labels      - {1 x M} x-tick labels              (default: '1','2',…)
%     colors      - [M x 3] RGB color per column        (default: lines(M))
%     jitter      - x-jitter half-width (data units)    (default: 0.15)
%     line_alpha  - opacity of connecting lines [0,1]   (default: 0.3)
%     marker_size - scatter dot area (pt^2)             (default: 40)
%     box         - underlay a neutral boxchart         (default: false)
%     stats       - pairwise paired tests with FDR      (default: false)
%     test        - 'signrank' | 'ttest'                (default: 'signrank')
%
%   Output:
%     h - struct with graphics handles:
%           .scatter  [M x 1] Scatter objects
%           .lines    {nLines x 1} cell of Line objects
%           .box      BoxChart array ([] when box=false)
%           .stats    annotation handles ([] when stats=false or n.s.)
%
%   See also: addSignificanceAnnotation, statsUtils, plotSplitViolin

arguments
    ax
    data   double
    opts.xpos        double  = []
    opts.labels      cell    = {}
    opts.colors              = []
    opts.jitter      double  = 0.15
    opts.line_alpha  double  = 0.3
    opts.marker_size double  = 40
    opts.box         logical = false
    opts.stats       logical = false
    opts.test        char {mustBeMember(opts.test, {'signrank','ttest'})} = 'signrank'
end

if isempty(ax); ax = gca; end

[V, M] = size(data);

if isempty(opts.xpos)
    opts.xpos = 1:M;
end
if isempty(opts.labels)
    opts.labels = arrayfun(@(m) num2str(m), 1:M, 'UniformOutput', false);
end
if isempty(opts.colors)
    opts.colors = lines(M);
end

hold_state = ishold(ax);
hold(ax, 'on');

h = struct('scatter', [], 'lines', {{}}, 'box', [], 'stats', []);

% --- Per-variable jitter: constant across columns so lines are clean diagonals ---
jitter_vals = (rand(V, 1) - 0.5) * 2 * opts.jitter;

% --- Connecting lines (drawn first, behind scatter) ---
for v = 1:V
    yv = data(v, :);
    if all(isnan(yv)); continue; end
    xv = opts.xpos + jitter_vals(v);
    hl = plot(ax, xv, yv, 'Color', [0.4 0.4 0.4], 'LineWidth', 0.7);
    hl.Color = [hl.Color(1:3), opts.line_alpha];
    h.lines{end+1} = hl;
end

% --- Scatter per column ---
h.scatter = gobjects(M, 1);
for m = 1:M
    xm = opts.xpos(m) + jitter_vals;
    ym = data(:, m);
    h.scatter(m) = scatter(ax, xm, ym, opts.marker_size, ...
        opts.colors(m, :), 'filled', ...
        'MarkerFaceAlpha', 0.85, 'MarkerEdgeColor', 'none');
end

% --- Overlay boxplot (RF_mkBoxPlot3 draws to current axes, no handles returned) ---
if opts.box
    datacells = num2cell(data, 1); % {1 x M}, each [V x 1]
    set(ancestor(ax, 'figure'), 'CurrentAxes', ax);
    RF_mkBoxPlot3(datacells, [0.5 0.5 0.5], opts.xpos, 0.4, 0, 1, 0, [0.88 0.88 0.88]);
end

% --- Axes cosmetics ---
ax.XTick      = opts.xpos;
ax.XTickLabel = opts.labels;
if M > 1
    margin = min(diff(sort(opts.xpos))) * 0.6;
else
    margin = 0.5;
end
xlim(ax, [min(opts.xpos) - margin, max(opts.xpos) + margin]);

% --- Paired stats ---
if opts.stats && M >= 2
    pairs  = nchoosek(1:M, 2);   % [nPairs x 2]
    nPairs = size(pairs, 1);
    p_raw  = nan(nPairs, 1);
    n_valid = zeros(nPairs, 1);

    for k = 1:nPairs
        d1 = data(:, pairs(k,1));
        d2 = data(:, pairs(k,2));
        ok = ~isnan(d1) & ~isnan(d2);
        n_valid(k) = sum(ok);
        if n_valid(k) < 3; continue; end
        switch opts.test
            case 'signrank'
                p_raw(k) = signrank(d1(ok), d2(ok));
            case 'ttest'
                [~, p_raw(k)] = ttest(d1(ok), d2(ok));
        end
    end

    p_adj = statsUtils.fdr(p_raw);

    % Build stats struct array (one entry per pair)
    stats(nPairs) = struct();
    for k = 1:nPairs
        stats(k).i           = pairs(k,1);
        stats(k).j           = pairs(k,2);
        stats(k).label_i     = opts.labels{pairs(k,1)};
        stats(k).label_j     = opts.labels{pairs(k,2)};
        stats(k).n           = n_valid(k);
        stats(k).test        = opts.test;
        stats(k).p_raw       = p_raw(k);
        stats(k).p_adj       = p_adj(k);
        stats(k).significant = p_adj(k) < 0.05;
    end
    h.stats = stats;

    % Annotate significant pairs; narrow spans sit lower
    xpos_pairs = opts.xpos(pairs);                             % [nPairs x 2] actual x
    [~, sort_idx] = sort(abs(xpos_pairs(:,2) - xpos_pairs(:,1)));

    yl           = ylim(ax);
    y_range      = yl(2) - yl(1);
    bracket_step = 0.08 * y_range;
    bracket_y    = yl(2) + bracket_step;
    any_annotated = false;

    for k = sort_idx'
        if isnan(p_adj(k)) || p_adj(k) >= 0.05; continue; end
        hs = addSignificanceAnnotation(ax, xpos_pairs(k,:), bracket_y, p_adj(k), ...
            'style', 'bracket');
        if ~isempty(hs)
            bracket_y     = bracket_y + bracket_step;
            any_annotated = true;
        end
    end

    if any_annotated
        ylim(ax, [yl(1), bracket_y + 0.5 * bracket_step]);
    end
end

if ~hold_state; hold(ax, 'off'); end

end
