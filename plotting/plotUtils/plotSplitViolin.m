function h = plotSplitViolin(ax, data1, data2, opts)
% plotSplitViolin  Split-violin plot of two distributions on given axes.
%
%   h = plotSplitViolin(ax, data1, data2)
%   h = plotSplitViolin(ax, data1, data2, Name, Value, ...)
%
%   Draws data1 as the left half and data2 as the right half of a split
%   violin centred at opts.x, with a small gap between the two halves.
%   Each half is normalised independently to its own peak (KDE or histogram),
%   so distribution shapes are directly comparable regardless of sample size.
%
%   Visualization mode:
%     - Default (edges=[])  : smooth KDE violin
%     - edges provided      : step-histogram bars using the supplied bin edges,
%                             applied identically to both distributions for
%                             direct visual comparison
%
%   Inputs:
%     ax    - axes handle, or [] to use gca
%     data1 - numeric vector, first distribution  (plotted on left)
%     data2 - numeric vector, second distribution (plotted on right)
%
%   Name-Value Options:
%     x            - x-coordinate of violin centre              (default: 1)
%     width        - half-width of each violin half (x units)   (default: 0.4)
%     gap          - gap between left and right halves (x units)(default: 0.04)
%     edges        - bin edges for histogram mode ([] = KDE)     (default: [])
%     color1       - fill colour for data1 / left half          (default: [0.20 0.47 0.76])
%     color2       - fill colour for data2 / right half         (default: [0.84 0.30 0.20])
%     fill_alpha   - patch face transparency [0,1]              (default: 0.7)
%     show_median  - draw median line                           (default: true)
%     show_mean    - draw mean line                             (default: false)
%     lw           - line width for median / mean lines         (default: 2)
%     median_color - colour of median lines                     (default: 'k')
%     mean_color   - colour of mean lines                       (default: 'r')
%     stats        - run kstest2 and annotate significance       (default: false)
%
%   Output:
%     h - struct of graphics handles:
%           .patch1, .patch2       filled violin / histogram patches
%           .median1, .median2     median lines  ([] when show_median=false)
%           .mean1, .mean2         mean lines    ([] when show_mean=false)
%           .stats                 annotation    ([] when stats=false or n.s.)
%
%   Examples:
%     % KDE mode
%     figure; ax = axes;
%     plotSplitViolin(ax, randn(50,1), randn(80,1)+1.5, 'stats', true);
%
%     % Histogram mode with common bin edges
%     edges = linspace(-4, 6, 20);
%     plotSplitViolin(ax, d1, d2, 'edges', edges, 'color1', [0.2 0.6 0.4]);
%
%   See also: addSignificanceAnnotation, statsUtils, ksdensity

    arguments
        ax
        data1 double {mustBeVector}
        data2 double {mustBeVector}
        opts.x           double  = 1
        opts.width       double  = 0.4
        opts.gap         double  = 0.04
        opts.edges               = []
        opts.color1               = [0.20 0.47 0.76]
        opts.color2               = [0.84 0.30 0.20]
        opts.fill_alpha  double  = 0.7
        opts.show_median logical = true
        opts.show_mean   logical = false
        opts.lw          double  = 2
        opts.median_color         = 'k'
        opts.mean_color           = 'r'
        opts.stats       logical = false
    end

    if isempty(ax)
        ax = gca;
    end

    data1 = data1(:);
    data2 = data2(:);

    % Preserve caller's hold state
    hold_state = ishold(ax);
    hold(ax, 'on');

    % Positions of the inner edges of each half
    xL = opts.x - opts.gap / 2;   % right boundary of left half
    xR = opts.x + opts.gap / 2;   % left  boundary of right half
    w  = opts.width;

    % ------------------------------------------------------------------ %
    % Build patches — KDE or histogram mode
    % ------------------------------------------------------------------ %
    if isempty(opts.edges)
        % ---- Smooth KDE mode ---------------------------------------- %
        [f1, y1] = ksdensity(data1);
        [f2, y2] = ksdensity(data2);

        % Independent normalisation: each peak reaches opts.width
        f1s = (f1 / max(f1)) * w;
        f2s = (f2 / max(f2)) * w;

        N1 = numel(y1);
        N2 = numel(y2);

        % Left patch: up the left edge, back down the inner boundary
        px1 = [xL - f1s,         repmat(xL, 1, N1)];
        py1 = [y1,                fliplr(y1)        ];
        h.patch1 = patch(ax, px1, py1, opts.color1, ...
            'FaceAlpha', opts.fill_alpha, 'EdgeColor', 'none');

        % Right patch: up the inner boundary, back down the right edge
        px2 = [repmat(xR, 1, N2), xR + fliplr(f2s)];
        py2 = [y2,                fliplr(y2)        ];
        h.patch2 = patch(ax, px2, py2, opts.color2, ...
            'FaceAlpha', opts.fill_alpha, 'EdgeColor', 'none');

        % References for median/mean interpolation and stats y_top
        y1_ref  = y1;    f1s_ref = f1s;
        y2_ref  = y2;    f2s_ref = f2s;
        y_top   = max([max(y1), max(y2)]);

    else
        % ---- Step-histogram mode ------------------------------------- %
        edges = opts.edges(:)';
        nb    = numel(edges) - 1;

        counts1 = histcounts(data1, edges);
        counts2 = histcounts(data2, edges);

        % Independent normalisation
        mx1 = max(counts1);  if mx1 == 0; mx1 = 1; end
        mx2 = max(counts2);  if mx2 == 0; mx2 = 1; end
        norm1 = (counts1 / mx1) * w;
        norm2 = (counts2 / mx2) * w;

        % Staircase y-coordinates shared by both sides:
        %   [e1, e2, e2, e3, e3, ..., eN]  (length 2*nb)
        inner  = reshape([edges(2:end-1); edges(2:end-1)], 1, []);
        y_step = [edges(1), inner, edges(end)];

        % Paired x-extents: [n1, n1, n2, n2, ..., nN]  (length 2*nb)
        x_n1 = reshape([norm1; norm1], 1, []);
        x_n2 = reshape([norm2; norm2], 1, []);

        % Left patch
        px1 = [xL - x_n1,         repmat(xL, 1, 2*nb)];
        py1 = [y_step,             fliplr(y_step)      ];
        h.patch1 = patch(ax, px1, py1, opts.color1, ...
            'FaceAlpha', opts.fill_alpha, 'EdgeColor', 'none');

        % Right patch
        px2 = [repmat(xR, 1, 2*nb), xR + fliplr(x_n2)];
        py2 = [y_step,              fliplr(y_step)     ];
        h.patch2 = patch(ax, px2, py2, opts.color2, ...
            'FaceAlpha', opts.fill_alpha, 'EdgeColor', 'none');

        % References for median/mean interpolation (bin centres)
        bin_centers = (edges(1:end-1) + edges(2:end)) / 2;
        y1_ref  = bin_centers;  f1s_ref = norm1;
        y2_ref  = bin_centers;  f2s_ref = norm2;
        y_top   = edges(end);
    end

    % ------------------------------------------------------------------ %
    % Median lines (from inner edge to violin/bar boundary at that y)
    % ------------------------------------------------------------------ %
    h.median1 = [];
    h.median2 = [];
    if opts.show_median
        m1 = median(data1, 'omitmissing');
        m2 = median(data2, 'omitmissing');

        w1_m = interp1(y1_ref, f1s_ref, m1, 'linear', 0);
        w2_m = interp1(y2_ref, f2s_ref, m2, 'linear', 0);

        h.median1 = plot(ax, [xL - w1_m, xL], [m1, m1], ...
            'Color', opts.median_color, 'LineWidth', opts.lw);
        h.median2 = plot(ax, [xR, xR + w2_m], [m2, m2], ...
            'Color', opts.median_color, 'LineWidth', opts.lw);
    end

    % ------------------------------------------------------------------ %
    % Mean lines
    % ------------------------------------------------------------------ %
    h.mean1 = [];
    h.mean2 = [];
    if opts.show_mean
        mu1 = mean(data1, 'omitmissing');
        mu2 = mean(data2, 'omitmissing');

        w1_mu = interp1(y1_ref, f1s_ref, mu1, 'linear', 0);
        w2_mu = interp1(y2_ref, f2s_ref, mu2, 'linear', 0);

        h.mean1 = plot(ax, [xL - w1_mu, xL], [mu1, mu1], ...
            'Color', opts.mean_color, 'LineWidth', opts.lw);
        h.mean2 = plot(ax, [xR, xR + w2_mu], [mu2, mu2], ...
            'Color', opts.mean_color, 'LineWidth', opts.lw);
    end

    % ------------------------------------------------------------------ %
    % KS test + significance bracket above the violin
    % ------------------------------------------------------------------ %
    h.stats = [];
    if opts.stats
        [~, p_ks] = kstest2(data1, data2);

        yl    = ylim(ax);
        y_ann = y_top + 0.04 * (yl(2) - yl(1));

        h.stats = addSignificanceAnnotation(ax, [xL - w, xR + w], ...
            y_ann, p_ks, 'style', 'bracket');
    end

    % Restore hold state
    if ~hold_state
        hold(ax, 'off');
    end
end
