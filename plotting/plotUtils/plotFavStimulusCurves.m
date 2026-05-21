function ax = plotFavStimulusCurves(quant, cfg, ax, normalize)
% PLOTFAVSTIMULUŚCURVES  Plot P(same favourite stimulus) vs anatomical distance.
%
%   ax = plotFavStimulusCurves(quant, cfg)
%   ax = plotFavStimulusCurves(quant, cfg, ax)
%   ax = plotFavStimulusCurves(quant, cfg, ax, true)   % normalized
%
%   quant      struct array [nSubjects x 1] from TopographicAnalysis.plotFavouriteStimulus
%   cfg        PlotConfig
%   ax         target axes; created in a new figure if omitted
%   normalize  logical (default false)
%              false : raw P(same fav), each curve's null shown as a dashed line
%                      of the same (lightened) colour
%              true  : z-score each curve against its own null — removes odour-
%                      dominance offsets so curves are comparable; null → y = 0
%
%   When numel(quant) > 1, subjects are pooled onto a common distance axis
%   (built from the union of all bin grids at the shared binSize).  Each
%   curve shows mean ± SEM across subjects; bins present in only a subset of
%   subjects contribute only those subjects to the mean.
    arguments
        quant     struct
        cfg       PlotConfig = PlotConfig()
        ax        = []
        normalize (1,1) logical = false
    end

    if isempty(ax)
        figure('Color', cfg.bgcol);
        ax = gca;
    end
    axes(ax);

    nSubjects = numel(quant);
    nReps     = size(quant(1).rep_curves, 2);
    binSize   = quant(1).binSize;

    % ---- Build common distance axis (union of all subjects' bin grids) ----
    all_min = arrayfun(@(q) min(q.binCenters), quant(:)');
    all_max = arrayfun(@(q) max(q.binCenters), quant(:)');
    c_min   = min(all_min);
    c_max   = max(all_max);
    x       = c_min : binSize : c_max;   % common bin-centre axis
    nCommon = numel(x);

    % Map each subject's curves to the common grid (exact, no interpolation
    % needed because bin centres are snapped to multiples of binSize).
    function mat = toCommon(vals, centers)
        % vals    [nBins_subj x nCols]
        % centers [1 x nBins_subj]
        nCols = size(vals, 2);
        mat   = nan(nCommon, nCols);
        % Find where this subject's bins sit in the common axis
        [~, ia] = ismembertol(centers(:), x(:), 1e-9);
        valid = ia > 0;
        mat(ia(valid), :) = vals(valid, :);
    end

    % ---- Pool curves onto common grid: [nCommon x nSubjects] per curve ----
    pool_rep  = @(r) cell2mat(arrayfun(@(q) toCommon(q.rep_curves(:,r),    q.binCenters), quant(:)', 'UniformOutput', false));
    pool_smr  = @(r) cell2mat(arrayfun(@(q) toCommon(q.rep_shuf_mean(:,r), q.binCenters), quant(:)', 'UniformOutput', false));
    pool_sst  = @(r) cell2mat(arrayfun(@(q) toCommon(q.rep_shuf_std(:,r),  q.binCenters), quant(:)', 'UniformOutput', false));
    pool_avg  =      cell2mat(arrayfun(@(q) toCommon(q.avg_curve,           q.binCenters), quant(:)', 'UniformOutput', false));
    pool_asm  =      cell2mat(arrayfun(@(q) toCommon(q.avg_shuf_mean,       q.binCenters), quant(:)', 'UniformOutput', false));
    pool_ass  =      cell2mat(arrayfun(@(q) toCommon(q.avg_shuf_std,        q.binCenters), quant(:)', 'UniformOutput', false));

    rep_colors = cfg.c(1:nReps, :);
    leg_h      = gobjects(nReps + 1, 1);
    lcfg_data  = struct('lineStyle', '-',  'lineWidth', 1.5);
    lcfg_avg   = struct('lineStyle', '--', 'lineWidth', 2.0);
    lcfg_null  = struct('lineStyle', '--', 'lineWidth', 0.8);

    if normalize
        % ---- Normalized: z-score per subject per bin, then pool ----
        % Null reference: y = 0 ± 1 σ band
        fill([x, fliplr(x)], [ones(1,nCommon), -ones(1,nCommon)], [0.6 0.6 0.6], ...
            'FaceAlpha', 0.15, 'EdgeColor', 'none');
        hold on;
        yline(0, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);

        for r = 1:nReps
            mat_r  = pool_rep(r);   shuf_mu_r = pool_smr(r);   shuf_sd_r = pool_sst(r);
            z_r    = (mat_r - shuf_mu_r) ./ max(shuf_sd_r, eps);   % [nCommon x nSubj]
            mu_z   = mean(z_r,  2, 'omitmissing');
            err_z  = std(z_r,   [], 2, 'omitmissing') / sqrt(nSubjects);
            leg_h(r) = plotLineNShade(x, mu_z, err_z, rep_colors(r,:), lcfg_data);
        end

        z_avg    = (pool_avg - pool_asm) ./ max(pool_ass, eps);
        mu_zavg  = mean(z_avg, 2, 'omitmissing');
        err_zavg = std(z_avg,  [], 2, 'omitmissing') / sqrt(nSubjects);
        leg_h(nReps + 1) = plotLineNShade(x, mu_zavg, err_zavg, cfg.axcol, lcfg_avg);

        ylabel(ax, 'P(same fav.) – null (σ)');

    else
        % ---- Raw: per-curve null as dashed lightened reference ----
        hold(ax, 'on');

        % Null dashed lines first (drawn beneath ribbons)
        for r = 1:nReps
            mu_null_r = mean(pool_smr(r), 2, 'omitmissing');
            null_col  = rep_colors(r,:) * 0.5 + 0.5;
            plot(ax, x, mu_null_r, '--', 'Color', null_col, 'LineWidth', lcfg_null.lineWidth);
        end
        mu_null_avg = mean(pool_asm, 2, 'omitmissing');
        plot(ax, x, mu_null_avg, '--', 'Color', [0.65 0.65 0.65], 'LineWidth', lcfg_null.lineWidth);

        % Real curves with ± SEM ribbons on top
        for r = 1:nReps
            mat_r  = pool_rep(r);
            mu_r   = mean(mat_r, 2, 'omitmissing');
            err_r  = std(mat_r,  [], 2, 'omitmissing') / sqrt(nSubjects);
            leg_h(r) = plotLineNShade(x, mu_r, err_r, rep_colors(r,:), lcfg_data);
        end

        mu_avg  = mean(pool_avg, 2, 'omitmissing');
        err_avg = std(pool_avg,  [], 2, 'omitmissing') / sqrt(nSubjects);
        leg_h(nReps + 1) = plotLineNShade(x, mu_avg, err_avg, cfg.axcol, lcfg_avg);

        ylabel(ax, 'P(same fav. stimulus)');
    end

    hold(ax, 'off');
    xlabel(ax, sprintf('Distance (%s)', distUnit(quant(1))));
    rep_leg = arrayfun(@(r) sprintf('Rep %d', r), 1:nReps, 'UniformOutput', false);
    legend(ax, leg_h, [rep_leg, {'Rep avg.'}], 'Location', 'best');
    set(ax, 'Color', cfg.bgcol, 'XColor', cfg.axcol, 'YColor', cfg.axcol);
    box(ax, 'off');
end

function u = distUnit(q)
% Heuristic: if max distance is very small it's probably µm-scale aligned data
    if max(q.binCenters) < 5000
        u = 'µm';
    else
        u = 'px';
    end
end
