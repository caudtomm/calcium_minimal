classdef VisualizeCellAnatomy
% VisualizeCellAnatomy  Visualize per-trial anatomy for top/bottom scored cells.
%
%   Picks n_cells from the top and n_cells from the bottom of a pre-computed
%   score distribution and plots, for each cell, its cropped anatomy image
%   with ROI boundary overlay across all filter-passing trials.
%   If scores are omitted the cells are picked at random.
%
%   Three figures are produced: anatomy_imgs, localcorr_imgs, and score summary.
%
%   Usage:
%       vca = VisualizeCellAnatomy(v);
%       vca.n_cells = 5;            % optional — default 5
%       vca.neighborhood = 25;      % optional — default 20 px
%       [hf_a, hf_lc, hf_dff, hf_s, anat, lcorr, dff, tbl] = vca.plot(scores);
%       [hf_a, hf_lc, hf_dff, hf_s, anat, lcorr, dff, tbl] = vca.plot();   % random cells
%
%   scores must be the {n_filtered_subjects} cell of [N_cells x 1] doubles
%   returned by BaselineDriftAnalysis.getContributionScores, using the same
%   ExperimentViewer v with the same DataFilter settings.
%
%   Outputs:
%       hf_anat    figure handle — anatomy images
%       hf_lcorr   figure handle — local correlation maps
%       hf_dff     figure handle — trial-averaged dF/F images (v.dataFilter.interval)
%       hf_scores  figure handle — regression lines + score distribution
%       anat_stacks  {n_rows} cell of [H_c x W_c x T_filtered]
%       lcorr_stacks {n_rows} cell of [H_c x W_c x T_filtered] ([] per row if unavailable)
%       dff_stacks   {n_rows} cell of [H_c x W_c x T_filtered] ([] per row if unavailable)
%       tbl        table with variables:
%                    row           figure row number
%                    subject_name  subject identifier string
%                    cell_id       1-based index into goodNeuron_IDs
%                    roi_id        ROI label in ROImap (NaN if file missing)
%                    score         drift score (NaN when scores not provided)

    properties
        v            ExperimentViewer
        n_cells      double = 5    % cells drawn from each end of score distribution
        neighborhood double = 20   % pixel padding around ROI (passed to ExtractCellAnatomy)
    end

    methods

        function obj = VisualizeCellAnatomy(v)
            arguments
                v ExperimentViewer
            end
            obj.v = v;
        end

        % ------------------------------------------------------------------
        function [hf_anat, hf_lcorr, hf_dff, hf_scores, anat_stacks, lcorr_stacks, dff_stacks, tbl] = plot(obj, scores)
        % plot  Main entry point.
        %
        %   [hf_anat, hf_lcorr, hf_scores, anat_stacks, lcorr_stacks, tbl] = plot()
        %   [hf_anat, hf_lcorr, hf_scores, anat_stacks, lcorr_stacks, tbl] = plot(scores)
            arguments
                obj    VisualizeCellAnatomy
                scores = {}
            end
            have_scores = ~isempty(scores) && iscell(scores);

            traces_list = obj.v.filtered_traces;   % {n_subj}
            n_subj      = numel(traces_list);

            %% -- Pool cells -----------------------------------------------
            all_s     = [];
            s_idx_vec = [];
            c_idx_vec = [];

            for i = 1:n_subj
                if have_scores && i <= numel(scores)
                    s = scores{i}(:);
                else
                    s = nan(traces_list{i}.N, 1);
                end
                n_cells_i  = numel(s);
                all_s      = [all_s;      s];                           %#ok<AGROW>
                s_idx_vec  = [s_idx_vec;  repmat(i, n_cells_i, 1)];    %#ok<AGROW>
                c_idx_vec  = [c_idx_vec;  (1:n_cells_i)'];              %#ok<AGROW>
            end

            % Keep only cells with finite scores (or all when no scores given)
            if have_scores
                valid = isfinite(all_s);
            else
                valid = true(size(all_s));
            end
            all_s     = all_s(valid);
            s_idx_vec = s_idx_vec(valid);
            c_idx_vec = c_idx_vec(valid);

            N_total = numel(all_s);
            n = min(obj.n_cells, floor(N_total / 2));
            if n == 0
                error('VisualizeCellAnatomy:notEnoughCells', ...
                    'Need at least 2 valid cells to visualise; found %d.', N_total);
            end

            %% -- Select cells ---------------------------------------------
            if have_scores
                [~, sort_ord] = sort(all_s);
                % top n: highest scores, shown in descending order (row 1 = highest)
                top_pos = flipud(sort_ord(end - n + 1 : end));
                % bottom n: lowest scores, shown in ascending order (row n+1 = most negative)
                bot_pos = sort_ord(1 : n);
            else
                perm    = randperm(N_total, 2 * n);
                top_pos = perm(1:n)';
                bot_pos = perm(n+1:end)';
            end
            sel_pos = [top_pos(:); bot_pos(:)];   % [2n x 1]
            n_rows  = numel(sel_pos);

            %% -- Filtered trial indices per subject -----------------------
            trial_idx_by_subj = cell(n_subj, 1);
            for i = 1:n_subj
                trial_idx_by_subj{i} = VisualizeCellAnatomy.getFilteredTrialIndices( ...
                    obj.v.dataFilter, traces_list{i});
            end

            % max_T: maximum over subjects that actually appear in sel_pos
            sel_subjs = unique(s_idx_vec(sel_pos));
            max_T = max(cellfun(@numel, trial_idx_by_subj(sel_subjs)));
            if max_T == 0
                error('VisualizeCellAnatomy:noTrials', ...
                    'No trials pass the current DataFilter settings.');
            end

            %% -- Extract anatomy and dF/F for each selected cell ----------
            subj_names = buildSubjectNames(obj.v, n_subj);

            anat_stacks   = cell(n_rows, 1);
            lcorr_stacks  = cell(n_rows, 1);
            dff_stacks    = cell(n_rows, 1);
            roi_masks     = cell(n_rows, 1);
            roi_masks_dff = cell(n_rows, 1);
            roi_ids       = nan(n_rows, 1);
            score_vals    = nan(n_rows, 1);
            interval_sec  = obj.v.dataFilter.interval;

            for r = 1:n_rows
                pos = sel_pos(r);
                si  = s_idx_vec(pos);
                ci  = c_idx_vec(pos);

                fprintf('\n[%d/%d] %s — cell %d\n', r, n_rows, subj_names{si}, ci);

                eca = ExtractCellAnatomy(traces_list{si}, ci, ...
                                         'neighborhood', obj.neighborhood);

                fprintf('  loading anatomy...');
                [anat_full, lcorr_full, roi_mask, roi_label] = eca.getCroppedStack();

                score_vals(r) = all_s(pos);

                %% Anatomy + localcorr ------------------------------------
                if ~isempty(anat_full)
                    fprintf(' %d trials\n', size(anat_full, 3));
                    t_idx = trial_idx_by_subj{si};
                    t_idx = t_idx(t_idx <= size(anat_full, 3));

                    anat_stacks{r}  = anat_full(:, :, t_idx);
                    roi_masks{r}    = roi_mask;
                    roi_ids(r)      = roi_label;

                    if ~isempty(lcorr_full)
                        t_idx_lc        = t_idx(t_idx <= size(lcorr_full, 3));
                        lcorr_stacks{r} = lcorr_full(:, :, t_idx_lc);
                    end
                else
                    fprintf(' unavailable\n');
                end

                %% dF/F mean images per trial -----------------------------
                try
                    [r1, r2, c1, c2, roi_mask_dff_r, roi_lbl_dff] = eca.getCropBounds();
                    t_idx_dff = trial_idx_by_subj{si};
                    n_t       = numel(t_idx_dff);
                    dff_imgs  = nan(r2-r1+1, c2-c1+1, n_t);
                    fprintf('  loading dF/F (%d trials): ', n_t);
                    for j = 1:n_t
                        fprintf('%d ', t_idx_dff(j));
                        fr_int = computeFrameInterval( ...
                            traces_list{si}, t_idx_dff(j), interval_sec);
                        snip = eca.loadTrialMovieSnippet(t_idx_dff(j), fr_int);
                        if ~isempty(snip) && ~isempty(snip.stack)
                            full_mean       = mean(snip.stack, 3, 'omitmissing');
                            dff_imgs(:,:,j) = full_mean(r1:r2, c1:c2);
                        end
                    end
                    fprintf('done\n');
                    dff_stacks{r}    = dff_imgs;
                    roi_masks_dff{r} = roi_mask_dff_r;
                    if isnan(roi_ids(r)); roi_ids(r) = roi_lbl_dff; end
                catch ME
                    fprintf('failed\n');
                    warning('VisualizeCellAnatomy:dffFailed', ...
                        'Could not extract dF/F images for row %d: %s', r, ME.message);
                end
            end

            %% -- Build figures --------------------------------------------
            fprintf('\nBuilding figures...\n');
            fprintf('  anatomy...'); drawnow;
            hf_anat   = obj.buildFigure(anat_stacks,  roi_masks,     max_T, 'Anatomy');
            fprintf(' done\n');
            fprintf('  local correlation...'); drawnow;
            hf_lcorr  = obj.buildFigure(lcorr_stacks, roi_masks,     max_T, 'Local correlation');
            fprintf(' done\n');
            fprintf('  dF/F mean...'); drawnow;
            hf_dff    = obj.buildFigure(dff_stacks,   roi_masks_dff, max_T, 'dF/F mean');
            fprintf(' done\n');
            fprintf('  score summary...'); drawnow;
            hf_scores = obj.buildScoreFigure(n, sel_pos, all_s, s_idx_vec, c_idx_vec, have_scores);
            fprintf(' done\n');

            %% -- Output table ---------------------------------------------
            row_nums  = (1:n_rows)';
            names_col = arrayfun(@(r) subj_names{s_idx_vec(sel_pos(r))}, ...
                                  1:n_rows, 'UniformOutput', false)';
            cell_ids  = arrayfun(@(r) c_idx_vec(sel_pos(r)), 1:n_rows)';

            tbl = table(row_nums, names_col, cell_ids, roi_ids, score_vals, ...
                'VariableNames', {'row', 'subject_name', 'cell_id', 'roi_id', 'score'});
        end

    end % public methods

    % ------------------------------------------------------------------
    methods (Access = private)

        function hf = buildFigure(obj, img_stacks, roi_masks, max_T, fig_title) %#ok<INUSL>
        % buildFigure  Render one cells-x-trials figure from a stack cell array.
        %
        %   img_stacks  {n_rows} of [H_c x W_c x T_r] or [] (blank row)
        %   roi_masks   {n_rows} of [H_c x W_c] or []
        %   max_T       number of columns
        %   fig_title   string displayed as figure name

            n_rows = numel(img_stacks);
            hf = figure('Color', 'k', 'Units', 'normalized', ...
                        'OuterPosition', [0 0 1 1], 'Name', fig_title);
            tiledlayout(n_rows, max_T, 'TileSpacing', 'none', 'Padding', 'tight');
            colormap(gray);

            for r = 1:n_rows
                imgs_row = img_stacks{r};

                if isempty(imgs_row)
                    % No data for this row — fill with black tiles
                    for j = 1:max_T
                        ax = nexttile;
                        set(ax, 'Color', 'k');
                        axis(ax, 'off');
                    end
                    continue
                end

                T_r      = size(imgs_row, 3);
                roi_mask = roi_masks{r};

                % Consistent clim across all trials for this cell
                finite_px = imgs_row(isfinite(imgs_row));
                if isempty(finite_px)
                    clim_range = [0, 1];
                else
                    clim_range = [min(finite_px), max(finite_px)];
                    if clim_range(1) >= clim_range(2)
                        clim_range(2) = clim_range(1) + eps;
                    end
                end

                % ROI boundary computed once per cell
                if ~isempty(roi_mask)
                    [boundaries, ~] = bwboundaries(roi_mask, 'noholes');
                else
                    boundaries = {};
                end

                for j = 1:max_T
                    ax = nexttile;
                    set(ax, 'Color', 'k');

                    if j <= T_r
                        imagesc(ax, imgs_row(:, :, j));
                        clim(ax, clim_range);
                        hold(ax, 'on');
                        for b = 1:numel(boundaries)
                            plot(ax, boundaries{b}(:, 2), boundaries{b}(:, 1), ...
                                 '-', 'Color', [1 1 0], 'LineWidth', 0.8);
                        end
                    end

                    axis(ax, 'image');
                    axis(ax, 'off');
                end
            end
        end

        % ------------------------------------------------------------------
        function hf = buildScoreFigure(obj, n, sel_pos, all_s, s_idx_vec, c_idx_vec, have_scores)
        % buildScoreFigure  Regression-line panel + score distribution histogram.
        %
        %   Layout: n rows x 3 columns.
        %     Col 1 : per-row baseline regression for high-score cells  (red)
        %     Col 2 : per-row baseline regression for low-score cells   (blue)
        %     Col 3 : pooled score histogram with xlines, spans all n rows

            col_hi = [0.80 0.15 0.10];
            col_lo = [0.10 0.30 0.80];

            B = obj.extractBaseline();   % {n_subj} of [N x T_trials]

            hf = figure('Color', 'w', 'Units', 'normalized', ...
                        'OuterPosition', [0 0 1 1], 'Name', 'Score summary');
            tiledlayout(n, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

            %% -- Col 3: score distribution (span all rows) ----------------
            ax_hist = nexttile(3, [n 1]);
            if have_scores
                histogram(ax_hist, all_s, 60, ...
                          'FaceColor', [0.72 0.72 0.72], 'EdgeAlpha', 0);
                hold(ax_hist, 'on');
                for r = 1:n
                    xline(ax_hist, all_s(sel_pos(r)), '-', 'Color', col_hi, ...
                          'LineWidth', 1.5, ...
                          'Label', sprintf('Row %d', r), ...
                          'LabelVerticalAlignment', 'bottom', 'FontSize', 7);
                    xline(ax_hist, all_s(sel_pos(n + r)), '-', 'Color', col_lo, ...
                          'LineWidth', 1.5, ...
                          'Label', sprintf('Row %d', n + r), ...
                          'LabelVerticalAlignment', 'bottom', 'FontSize', 7);
                end
                xlabel(ax_hist, 'Score (\sigma / trial)');
                ylabel(ax_hist, 'Cell count');
                title(ax_hist, 'Score distribution');
            else
                axis(ax_hist, 'off');
                text(0.5, 0.5, 'No scores provided', ...
                     'Units', 'normalized', 'HorizontalAlignment', 'center', ...
                     'Parent', ax_hist);
            end
            box(ax_hist, 'off');

            %% -- Cols 1 & 2: regression plots ----------------------------
            % nexttile without index fills tiles in reading order,
            % skipping the span region already occupied by ax_hist.
            for r = 1:n
                % High-score cell — anatomy figure row r
                pos_hi = sel_pos(r);
                si = s_idx_vec(pos_hi);  ci = c_idx_vec(pos_hi);
                ax_hi = nexttile;
                if si <= numel(B) && ~isempty(B{si}) && ci <= size(B{si}, 1)
                    plotRegression(ax_hi, B{si}(ci, :), col_hi, r, all_s(pos_hi));
                else
                    axis(ax_hi, 'off');
                end

                % Low-score cell — anatomy figure row n+r
                pos_lo = sel_pos(n + r);
                si = s_idx_vec(pos_lo);  ci = c_idx_vec(pos_lo);
                ax_lo = nexttile;
                if si <= numel(B) && ~isempty(B{si}) && ci <= size(B{si}, 1)
                    plotRegression(ax_lo, B{si}(ci, :), col_lo, n + r, all_s(pos_lo));
                else
                    axis(ax_lo, 'off');
                end
            end
        end

        % ------------------------------------------------------------------
        function B = extractBaseline(obj)
        % extractBaseline  Per-trial mean-baseline matrix for each subject.
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

    end % private methods

    methods (Access = private, Static)

        function trial_idx = getFilteredTrialIndices(df, traces)
        % Replicate DataFilter.filterData trial-index logic without loading
        % the actual trace data.  Returns chronological trial indices that
        % pass trial_nums, stim, and repetition filters.

            % Start chronological
            trial_idx = (1:traces.ntrials)';

            % --- trial_nums filter ---
            if ~isempty(df.trial_nums)
                trial_idx = trial_idx(ismember(trial_idx, df.trial_nums));
            end
            if isempty(trial_idx); return; end

            % --- stim filter ---
            stims = traces.stim_series.stimulus(trial_idx);
            try
                desired = getStimuliByGroup(traces.subject_group, df.stims_allowed);
            catch
                desired = {};   % unknown group / tag: keep all
            end
            if ~isempty(desired)
                keep      = ismember(stims, desired);
                trial_idx = trial_idx(keep);
                stims     = stims(keep);
            end
            if isempty(trial_idx); return; end

            % --- repetition filter ---
            reps_touse = df.repetitions;
            if isempty(reps_touse); return; end   % empty → keep all

            thisstims  = unique(stims);
            idx_keep   = false(numel(stims), 1);
            for s = 1:numel(thisstims)
                pos_s          = find(strcmp(stims, thisstims{s}));
                reps_available = 1:numel(pos_s);
                reps_valid     = reps_available(ismember(reps_available, reps_touse));
                if ~isempty(reps_valid)
                    idx_keep(pos_s(reps_valid)) = true;
                end
            end
            trial_idx = trial_idx(idx_keep);
        end

    end % private static methods
end % classdef


% ======================================================================
% File-private helper
% ======================================================================

function plotRegression(ax, y, color, row_num, score_val)
% plotRegression  Scatter + OLS fit for a single cell's per-trial baseline.
%
%   ax        axes handle
%   y         [1 x T] baseline values (may contain NaN)
%   color     [1x3] RGB line/marker colour
%   row_num   integer label for the anatomy figure row
%   score_val scalar drift score (may be NaN)
    T = numel(y);
    t = (1:T)';
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

    scatter(ax, t, y(:), 12, color, 'filled', 'MarkerFaceAlpha', 0.5);
    hold(ax, 'on');
    plot(ax, t, y_fit, '-', 'Color', color * 0.65, 'LineWidth', 1.5);
    xlabel(ax, 'Trial');
    ylabel(ax, 'Baseline (a.u.)');
    if isnan(score_val)
        title(ax, sprintf('Row %d', row_num), 'FontSize', 8);
    else
        title(ax, sprintf('Row %d   score = %.3f', row_num, score_val), 'FontSize', 8);
    end
    axis(ax, 'tight');
    box(ax, 'off');
end

% ======================================================================

function frame_interval = computeFrameInterval(traces, trial_num, interval_sec)
% computeFrameInterval  Convert a [t_start, t_end] interval (seconds relative
%   to effective stimulus onset) into movie frame indices.
%
%   Effective stimulus onset = stim_series.frame_onset + odor_delay * fs.
%   trial_num is used as a row index into stim_series (clamped to valid range).
    fs            = traces.framerate;
    row           = min(trial_num, height(traces.stim_series));
    stim_onset_fr = traces.stim_series.frame_onset(row) + round(traces.odor_delay * fs);
    fr_start      = stim_onset_fr + round(interval_sec(1) * fs);
    fr_end        = stim_onset_fr + round(interval_sec(2) * fs);
    frame_interval = fr_start : fr_end;
end

% ======================================================================

function names = buildSubjectNames(v, n_subj)
% Return {n_subj} cell of subject name strings for the filtered subjects.
    names = arrayfun(@(k) num2str(k), 1:n_subj, 'UniformOutput', false);

    if ~any(strcmp('name', v.subjectTab.Properties.VariableNames))
        % Fallback: use subject_ID from each filtered trace's Locations
        filtered_traces = v.filtered_traces;
        for k = 1:min(n_subj, numel(filtered_traces))
            sid = filtered_traces{k}.subject_locations.subject_ID;
            if ~isempty(sid)
                names{k} = sid;
            end
        end
        return
    end

    % Preferred: names from subjectTab, mapped through subjects_to_use
    global_idx = find(v.subjects_to_use);
    for k = 1:min(n_subj, numel(global_idx))
        nm = v.subjectTab.name(global_idx(k));
        if iscell(nm),  nm = nm{1}; end
        names{k} = char(nm);
    end
end
