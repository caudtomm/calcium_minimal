classdef BaselineDriftAnalysis
% BaselineDriftAnalysis  Score cells by contribution to baseline drift.
%
%   Wraps an ExperimentViewer and uses its current dataFilter to extract
%   per-trial baseline vectors, then scores each cell by how strongly it
%   drives the population-level drift across trials.
%
%   Usage:
%       bda    = BaselineDriftAnalysis(v);        % v: ExperimentViewer
%       scores = bda.getContributionScores('linear');
%       bda.plotScoreQC(scores, 'linear');
%
%   Methods:   'linear'        regression slope vs trial index
%              'general'       first PC of trial-to-trial displacements
%              'arch sum'      arch path with L1 (total-FR) constraint
%              'arch geodesic' arch path with L2 (hypersphere) constraint
%
%   Output:
%       scores  {n_subjects} cell of [N_cells x 1] signed, normalised
%               drift scores in units of [sigma / trial].

    properties
        v ExperimentViewer
    end

    methods (Access = public)

        function obj = BaselineDriftAnalysis(v)
            arguments
                v ExperimentViewer
            end
            obj.v = v;
        end

        function scores = getContributionScores(obj, method)
        % getContributionScores  Compute per-cell baseline drift scores.
        %
        %   scores = getContributionScores(method)
        %
        %   Input:
        %       method  char — one of 'linear', 'general',
        %                      'arch sum', 'arch geodesic'
        %
        %   Output:
        %       scores  {n_subjects} cell of [N_cells x 1] doubles.
        %               Sign: positive = cell activity increases in the
        %               direction of the dominant drift.

            B = obj.extractBaseline();   % {n_subj} of [N x T_trials]

            switch lower(strtrim(method))
                case 'linear'
                    scores = obj.scoreLinear(B);
                case 'general'
                    scores = obj.scoreGeneral(B);
                case 'arch sum'
                    scores = obj.scoreArch(B, 'sum');
                case 'arch geodesic'
                    scores = obj.scoreArch(B, 'geodesic');
                otherwise
                    error('BaselineDriftAnalysis:unknownMethod', ...
                        'Unknown method ''%s''. Choose: linear, general, arch sum, arch geodesic.', method);
            end
        end

        function hf = plotScoreQC(obj, scores, method)
        % plotScoreQC  QC figure for getContributionScores output.
        %
        %   hf = plotScoreQC(scores)
        %   hf = plotScoreQC(scores, method)   % method string for title
        %
        %   Layout (3 rows x 4 cols):
        %     Row 1 : three example-cell regressions + score histogram
        %     Row 2 : population heatmap (cells sorted by score x trials)
        %     Row 3 : mean drift trajectory (top/mid/bot deciles) +
        %             per-subject boxplot with pairwise Mann-Whitney

            narginchk(2, 3);
            if nargin < 3, method = ''; end

            %% -- Pool and sort scores ------------------------------------
            [all_s, s_idx, c_idx] = poolScores(scores);
            N_valid = numel(all_s);
            [~, sort_ord] = sort(all_s);

            %% -- Example cells: 98th / 50th / 2nd percentile -----------
            pct      = [0.98, 0.50, 0.02];
            ex_rank  = max(1, min(N_valid, round(pct * N_valid)));
            ex_glob  = sort_ord(ex_rank);   % indices into pooled valid array
            ex_labels = {'High +', 'Near median', 'High -'};
            ex_col    = {[0.80 0.15 0.10], [0.50 0.50 0.50], [0.10 0.30 0.80]};

            %% -- Baseline data ------------------------------------------
            B = obj.extractBaseline();      % {n_subj} of [N x T_trials]

            %% -- Figure -------------------------------------------------
            hf = figure('Color', 'w', 'Units', 'normalized', ...
                        'OuterPosition', [0.02 0.05 0.96 0.90]);
            tl = tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
            if ~isempty(method)
                title(tl, ['Baseline drift QC — ', method], 'FontWeight', 'bold');
            end

            %% -- Row 1: example cell regressions (tiles 1-3) -----------
            for k = 1:3
                ax = nexttile(k);
                g  = ex_glob(k);
                si = s_idx(g);  ci = c_idx(g);
                y  = double(B{si}(ci, :));   % [1 x T], may contain NaN  
                T  = numel(y);
                t  = (1:T)';

                % Use only valid trials — consistent with olsSlope                                                                                                                                                                 
                valid = ~isnan(y(:));                                                                                                                                                                                              
                if sum(valid) >= 2                                                                                                                                                                                                 
                   t_v   = t(valid);  y_v = y(valid);                                                                                                                                                                             
                   t_c_v = t_v - mean(t_v);                                                                                                                                                                                       
                   slope = (y_v(:)' * t_c_v) / sum(t_c_v .^ 2);                                                                                                                                                                   
                   intrc = mean(y_v) - slope * mean(t_v);                                                                                                                                                                         
                else                                                                                                                                                                                                               
                   slope = 0;                                                                                                                                                                                                     
                   intrc = mean(y, 'omitmissing');                                                                                                                                                                                
                end                                 
                y_fit = intrc + slope .* t;

                scatter(ax, t, y, 18, ex_col{k}, 'filled', ...
                        'MarkerFaceAlpha', 0.55);
                hold(ax, 'on');
                plot(ax, t, y_fit, '-', 'Color', ex_col{k} .* 0.65, ...
                     'LineWidth', 2);
                xlabel(ax, 'Trial');
                ylabel(ax, 'Baseline (a.u.)');
                title(ax, sprintf('%s  (subj %d, cell %d)\nscore = %.3f', ...
                    ex_labels{k}, si, ci, all_s(g)));
                axis(ax, 'tight');  box(ax, 'off');
            end

            %% -- Tile 4: score histogram --------------------------------
            ax4 = nexttile(4);
            histogram(ax4, all_s, 60, 'FaceColor', [0.72 0.72 0.72], ...
                      'EdgeAlpha', 0);
            hold(ax4, 'on');
            for k = 1:3
                xline(ax4, all_s(ex_glob(k)), '-', 'Color', ex_col{k}, ...
                      'LineWidth', 2, 'Label', ex_labels{k}, ...
                      'LabelVerticalAlignment', 'bottom', 'FontSize', 8);
            end
            xlabel(ax4, 'Score (\sigma / trial)');
            ylabel(ax4, 'Cell count');
            title(ax4, 'Score distribution');
            box(ax4, 'off');

            %% -- Row 2: population heatmap (tiles 5-8) -----------------
            ax_heat = nexttile(5, [1 4]);

            % Build [N_valid x max_T] z-scored matrix, rows sorted by score
            T_list = cellfun(@(b) size(b, 2), B);
            max_T  = max(T_list);
            Z = nan(N_valid, max_T);
            for k = 1:N_valid
                g  = sort_ord(k);
                si = s_idx(g);  ci = c_idx(g);
                y  = double(B{si}(ci, 1:T_list(si)));
                mu = nanmean(y);
                sg = nanstd(y);
                if sg < eps, continue; end
                Z(k, 1:T_list(si)) = (y - mu) ./ sg;
            end

            imagesc(ax_heat, 1:max_T, 1:N_valid, Z);
            colormap(ax_heat, obj.v.plotConfig.colormapName);
            finite_z = Z(isfinite(Z));
            c_rng    = max(eps, prctile(abs(finite_z), 95));   % guard empty / zero     
            clim(ax_heat, [-c_rng, c_rng]);
            cb = colorbar(ax_heat);
            cb.Label.String = 'z-score';
            xlabel(ax_heat, 'Trial');
            ylabel(ax_heat, 'Cell rank (low \rightarrow high score)');
            title(ax_heat, 'Population baseline activity');

            % Mark decile boundaries
            n10 = round(0.10 * N_valid);
            yline(ax_heat, n10,           '--', 'Color', ex_col{3}, ...
                  'LineWidth', 1, 'Alpha', 0.8);
            yline(ax_heat, N_valid - n10, '--', 'Color', ex_col{1}, ...
                  'LineWidth', 1, 'Alpha', 0.8);

            % Swap x and y axes
            view([-90 90])

            %% -- Row 3a: mean drift trajectory (tiles 9-11) ------------
            ax_traj = nexttile(9, [1 3]);

            n10      = max(1, round(0.10 * N_valid));
            top_rows = Z(end - n10 + 1 : end,  :);
            mid_rows = Z(max(1, round(N_valid*0.45)) : min(N_valid, round(N_valid*0.55)), :);
            bot_rows = Z(1 : n10,               :);

            top_mu = mean(top_rows, 1, 'omitmissing');
            mid_mu = mean(mid_rows, 1, 'omitmissing');
            bot_mu = mean(bot_rows, 1, 'omitmissing');

            t_ax = 1:max_T;
            plot(ax_traj, t_ax, top_mu, '-', 'Color', ex_col{1}, 'LineWidth', 2);
            hold(ax_traj, 'on');
            plot(ax_traj, t_ax, mid_mu, '-', 'Color', ex_col{2}, 'LineWidth', 1.5);
            plot(ax_traj, t_ax, bot_mu, '-', 'Color', ex_col{3}, 'LineWidth', 2);
            yline(ax_traj, 0, ':', 'Color', [0.6 0.6 0.6]);
            legend(ax_traj, {'Top 10%', 'Mid 10%', 'Bot 10%'}, ...
                   'Location', 'best', 'Box', 'off');
            xlabel(ax_traj, 'Trial');
            ylabel(ax_traj, 'Mean z-score');
            title(ax_traj, 'Mean drift trajectory by score decile');
            axis(ax_traj, 'tight');  box(ax_traj, 'off');

            %% -- Row 3b: per-subject comparison (tile 12) --------------
            ax_subj = nexttile(12);

            % Keep only subjects with >= 5 valid scores
            keep      = cellfun(@(s) sum(~isnan(s(:))) >= 5, scores);
            subj_data = cellfun(@(s) s(~isnan(s(:))), ...
                            scores(keep), 'UniformOutput', false);
            n_keep    = numel(subj_data);

            if n_keep == 0
                title(ax_subj, 'No subjects with >= 5 cells');
                return
            end

            axes(ax_subj); %#ok<LAXES>  RF_mkBoxPlot3 targets current axes
            RF_mkBoxPlot3(subj_data, [], [], .3, .5, .5, 4, []);
            ax_subj = gca;

            xlabel(ax_subj, 'Subject');
            ylabel(ax_subj, 'Score (\sigma / trial)');
            title(ax_subj, 'Per-subject distribution');
            xticks(ax_subj, 1:n_keep);

            % Subject labels
            subj_names = arrayfun(@(i) num2str(i), find(keep), ...
                                  'UniformOutput', false);
            if any(strcmp('name', obj.v.subjectTab.Properties.VariableNames)) && ~isempty(obj.v.subjectTab.name)
                subj_names = obj.v.subjectTab.name(keep);
            end
            % xticklabels(ax_subj, subj_names); % annoyingly long names
            % xtickangle(ax_subj, 45);
            box(ax_subj, 'off');

            % Pairwise Mann-Whitney with FDR
            if n_keep >= 2
                mw  = statsUtils.pairwiseMannWhitney(subj_data);
                sig = mw([mw.significant]);

                if ~isempty(sig)
                    % Sort by p_adj, annotate at most 8 brackets
                    [~, p_ord] = sort([sig.p_adj]);
                    sig        = sig(p_ord(1 : min(8, end)));

                    y_range   = diff(ylim(ax_subj));
                    bracket_y = max(ylim(ax_subj)) + 0.05 * y_range;
                    b_step    = 0.08 * y_range;

                    for m = 1:numel(sig)
                        addSignificanceAnnotation(ax_subj, ...
                            [sig(m).i, sig(m).j], bracket_y, ...
                            sig(m).p_adj, 'style', 'bracket');
                        bracket_y = bracket_y + b_step;
                    end
                end
            end
        end

    end % public methods

    % ------------------------------------------------------------------
    methods (Access = private)

        function B = extractBaseline(obj)
        % Extract mean-over-time baseline vector per trial.
        % Forces chronological trial ordering; restores original after.
        % Returns {n_subj} cell of [N_cells x N_trials].

            prev_sorting = obj.v.dataFilter.trial_sorting;
            obj.v.dataFilter.trial_sorting = 'chronological';

            [~, events] = ModeSelector(obj.v).extract;
            % events{i}: [T x N x trials]  ->  mean over T  ->  [N x trials]
            B = cellfun( ...
                @(x) reshape(mean(x, 1, 'omitmissing'), size(x,2), size(x,3)), ...
                events, 'UniformOutput', false);

            obj.v.dataFilter.trial_sorting = prev_sorting;
        end

        % --------------------------------------------------------------
        function scores = scoreLinear(~, B)
        % 'linear': per-cell OLS slope vs trial index, normalised by std.
        %
        % Score = slope_n / sigma_n  [sigma / trial]
        % Sign   = direction of monotonic trend.

            n_subj = numel(B);
            scores = cell(n_subj, 1);
            for i = 1:n_subj
                scores{i} = computeLinearScore(B{i});
            end
        end

        % --------------------------------------------------------------
        function scores = scoreGeneral(~, B)
        % 'general': first PC of trial-to-trial displacement vectors.
        %
        % Steps:
        %   1. Compute delta_t = b_{t+1} - b_t  [N x (T-1)]
        %   2. SVD of delta to find dominant displacement direction u [N x 1]
        %   3. mean_proj = mean projection of delta columns onto u  [scalar]
        %   4. score_n = u_n * mean_proj / sigma_n  [sigma / trial]
        %
        % Sign convention: u is flipped so that mean_proj >= 0
        % (i.e., the first PC points in the direction of average drift).
        %
        % Missing values: NaN trials are mean-imputed per cell before SVD.
        % Cells with < 2 valid trials are zeroed out (no drift contribution)
        % and receive NaN scores via sigma = NaN.

            n_subj = numel(B);
            scores = cell(n_subj, 1);
            for i = 1:n_subj
                b     = B{i};          % [N x T]
                if size(b, 2) < 2
                    scores{i} = nan(size(b, 1), 1);
                    continue;
                end
                sigma = std(b, [], 2, 'omitmissing');
                sigma = safeStd(sigma);

                % Impute NaN with per-cell mean so diff/SVD are well-defined.
                % Cells with < 2 valid trials are zeroed (contribute no drift).
                b_imp = imputeMean(b);
                delta = diff(b_imp, 1, 2); % [N x (T-1)]

                % First left singular vector = dominant displacement direction
                [U, ~, ~] = svd(delta, 'econ');
                u = U(:, 1);           % [N x 1], unit norm

                % Align sign with average drift direction
                mean_delta = mean(delta, 2);   % [N x 1], no NaN after imputation
                if dot(u, mean_delta) < 0
                    u = -u;
                end

                % Average per-trial displacement along drift direction
                mean_proj = mean(u' * delta); % [activity / trial]

                scores{i} = u .* mean_proj ./ sigma;
            end
        end

        % --------------------------------------------------------------
        function scores = scoreArch(~, B, archType)
        % 'arch sum' / 'arch geodesic':
        %   Derive linear-regression endpoints A (t=1) and Bend (t=T),
        %   interpolate a constrained path between them, then score each
        %   cell by the OLS slope of its arch-path coordinate vs trial.
        %
        %   arch sum:      straight line on L1-normalised (proportional)
        %                  vectors, scaled by linearly-interpolated total FR.
        %   arch geodesic: SLERP on L2 sphere, with linearly-interpolated norm.
        %
        %   Score = arch_slope_n / sigma_n_observed  [sigma / trial]
        %
        %   Missing values: path is built only from cells with valid endpoints
        %   (>= 2 non-NaN trials). Invalid cells receive NaN scores.

            n_subj = numel(B);
            scores = cell(n_subj, 1);
            for i = 1:n_subj
                b = B{i};                      % [N x T]
                N = size(b, 1);
                T = size(b, 2);
                sigma = std(b, [], 2, 'omitmissing');
                sigma = safeStd(sigma);

                [A_all, Bend_all] = linearEndpoints(b);  % NaN for cells with < 2 valid

                % Build path only from cells with valid (finite) endpoints.
                % A single NaN endpoint would corrupt norm/dot-product computations.
                valid = isfinite(A_all);       % [N x 1]
                arch_slope = nan(N, 1);

                if any(valid)
                    switch archType
                        case 'geodesic'
                            path_v = computeSlerpPath(A_all(valid), Bend_all(valid), T);
                        case 'sum'
                            path_v = computeSumPath(A_all(valid), Bend_all(valid), T);
                    end
                    arch_slope(valid) = olsSlope(path_v);  % path is complete, no NaN
                end

                scores{i} = arch_slope ./ sigma;
            end
        end

    end % private methods
end % classdef


% ======================================================================
% File-private helpers
% ======================================================================

function score = computeLinearScore(b)
% OLS slope of each cell's baseline activity vs trial index.
% b: [N x T].  Returns [N x 1] signed scores [sigma / trial].
    sigma = std(b, [], 2, 'omitmissing');
    sigma = safeStd(sigma);
    slope = olsSlope(b);       % [N x 1]
    score = slope ./ sigma;
end

% ----------------------------------------------------------------------
function slope = olsSlope(b)
% Vectorised NaN-aware OLS slope of rows of b against trial index 1..T.
% b: [N x T].  Returns slope [N x 1].
% Each row uses only its own non-NaN columns; rows with < 2 valid trials
% return NaN.
    T      = size(b, 2);
    t      = (1:T)';                          % [T x 1]
    mask   = ~isnan(b);                       % [N x T]
    T_n    = sum(mask, 2);                    % [N x 1] valid trial count

    b_safe = b;  b_safe(~mask) = 0;           % replace NaN with 0 for dot products
    t_mean = (mask * t)              ./ T_n;  % [N x 1] per-cell mean of valid t
    b_mean = sum(b_safe .* mask, 2)  ./ T_n;  % [N x 1] per-cell mean of valid b

    t_c  = (t' - t_mean) .* mask;             % [N x T] centred t, invalid zeroed
    b_c  = (b_safe - b_mean) .* mask;         % [N x T] centred b, invalid zeroed

    slope = sum(b_c .* t_c, 2) ./ sum(t_c .^ 2, 2);  % [N x 1]
    slope(T_n < 2) = NaN;
end

% ----------------------------------------------------------------------
function [A, Bend] = linearEndpoints(b)
% Predicted baseline vectors at t=1 and t=T from NaN-aware OLS regression.
% b: [N x T].  Returns A, Bend: [N x 1].  NaN for cells with < 2 valid trials.
    T      = size(b, 2);
    t      = (1:T)';
    mask   = ~isnan(b);
    T_n    = sum(mask, 2);
    b_safe = b;  b_safe(~mask) = 0;
    t_mean = (mask * t)              ./ T_n;   % [N x 1] per-cell mean of valid t
    b_mean = sum(b_safe .* mask, 2)  ./ T_n;   % [N x 1]
    slope  = olsSlope(b);                      % [N x 1], NaN for T_n < 2
    intercept = b_mean - slope .* t_mean;
    A    = intercept + slope .* 1;
    Bend = intercept + slope .* T;
end

% ----------------------------------------------------------------------
function path = computeSlerpPath(A, Bend, T)
% SLERP from A to Bend in T steps.
% Direction interpolated on unit sphere; norm interpolated linearly.
% Returns [N x T] path matrix.

    nA = norm(A);
    nB = norm(Bend);

    if nA < eps || nB < eps
        % Degenerate: fall back to linear
        path = A + (Bend - A) .* linspace(0, 1, T);
        return
    end

    ahat = A    / nA;
    bhat = Bend / nB;

    cosTheta = max(-1, min(1, dot(ahat, bhat)));
    theta    = acos(cosTheta);

    s     = linspace(0, 1, T);              % [1 x T]
    norms = nA + s .* (nB - nA);           % linearly interpolated norms

    if abs(sin(theta)) < 1e-10
        % A and Bend nearly parallel: linear direction interpolation
        dirs = ahat + s .* (bhat - ahat);  % [N x T]
        dirs = dirs ./ vecnorm(dirs, 2, 1);
    else
        dirs = (sin((1-s).*theta) .* ahat + sin(s.*theta) .* bhat) ...
               ./ sin(theta);              % [N x T]
    end

    path = dirs .* norms;                  % [N x T]
end

% ----------------------------------------------------------------------
function path = computeSumPath(A, Bend, T)
% L1-constrained path from A to Bend in T steps.
% Proportional direction straight-line on normalised vectors;
% total firing rate linearly interpolated.
% Returns [N x T] path matrix.

    frA = sum(A);
    frB = sum(Bend);

    if abs(frA) < eps || abs(frB) < eps
        warning('BaselineDriftAnalysis:zeroTotalFR', ...
            'Near-zero total FR at a regression endpoint; falling back to linear path.');
        path = A + (Bend - A) .* linspace(0, 1, T);
        return
    end

    ahat = A    / frA;
    bhat = Bend / frB;

    s        = linspace(0, 1, T);              % [1 x T]
    frs      = frA + s .* (frB - frA);        % linearly interpolated total FR
    dirs     = ahat + s .* (bhat - ahat);     % [N x T], straight line on proportions
    dir_sums = sum(dirs, 1);                  % [1 x T]

    % guard against degenerate columns
    bad = abs(dir_sums) < eps;
    dir_sums(bad) = 1;

    path = (dirs ./ dir_sums) .* frs;    % scale proportions back to activity
end

% ----------------------------------------------------------------------
function b_imp = imputeMean(b)
% Replace NaN entries with the per-row mean.
% Rows with < 2 valid (non-NaN) trials are zeroed out entirely so they
% contribute no drift signal to downstream SVD computations.
    T_n   = sum(~isnan(b), 2);              % [N x 1]
    b_imp = b;

    % Zero out rows that have too few valid trials
    b_imp(T_n < 2, :) = 0;

    % Fill remaining NaN positions with the row mean
    b_mean    = mean(b, 2, 'omitmissing');           % [N x 1]
    fill      = isnan(b_imp) & (T_n >= 2);           % [N x T]
    if any(fill(:))
        fill_vals = repmat(b_mean, 1, size(b, 2));   % [N x T]
        b_imp(fill) = fill_vals(fill);
    end
end

% ----------------------------------------------------------------------
function sigma = safeStd(sigma)
% Replace near-zero std with NaN to avoid division blow-up.
    sigma(sigma < eps) = NaN;
end

% ----------------------------------------------------------------------
function [all_s, s_idx, c_idx] = poolScores(scores)
% Pool {n_subj} score cells into flat valid vectors with origin tracking.
%
%   all_s  [N_valid x 1]  finite scores
%   s_idx  [N_valid x 1]  subject index for each score
%   c_idx  [N_valid x 1]  cell index within that subject
    all_s = []; s_idx = []; c_idx = [];
    for i = 1:numel(scores)
        s = scores{i}(:);
        n = numel(s);
        all_s = [all_s; s];                       %#ok<AGROW>
        s_idx = [s_idx; repmat(i, n, 1)];         %#ok<AGROW>
        c_idx = [c_idx; (1:n)'];                  %#ok<AGROW>
    end
    valid = isfinite(all_s);
    all_s = all_s(valid);
    s_idx = s_idx(valid);
    c_idx = c_idx(valid);
end
