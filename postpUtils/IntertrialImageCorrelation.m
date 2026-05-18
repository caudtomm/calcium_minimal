classdef IntertrialImageCorrelation
% IntertrialImageCorrelation  Intertrial image correlation for dF/F and anatomy.
%
% For each filtered subject, loads the trial-averaged dF/F image stack
% (from data_dir) and the anatomical image stack (from fish1.mat), crops
% each around every ROI, and computes intertrial Pearson correlation
% matrices.  The sheared correlation profiles for both modalities are
% plotted together as line + shade curves.
%
% Cells can be optionally filtered by their baseline drift score before
% the correlation is computed.  Scores are accepted as {n_subjects} cell
% of [N x 1] doubles (output of BaselineDriftAnalysis.getContributionScores)
% and the included population is defined by a quantile range on the pooled
% score distribution.
%
% Pipeline:
%   (1) Load per-subject [H x W x nTrials] stacks:
%         dF/F    — from data_dir, matched by subject_ID
%         anatomy — from subject file (subject.anatomy_imgs)
%   (2) For each ROI passing the score filter, crop both stacks with
%       neighborhood-pixel padding around the ROI bounding box.
%   (3) Compute the [nTrials x nTrials] Pearson intertrial correlation
%       matrix from the flattened crop (NaN-safe, 'complete' rows).
%   (4) Concatenate all subjects/ROIs along the 3rd dim →
%       [nTrials x nTrials x N_total]  (one per modality).
%   (5) Shear: shift column i up by (i-1), embedding into
%       [2*nTrials-1 x nTrials x N_total].  Diagonal lands at row nTrials.
%   (6) Plot both modalities as mean ± SEM vs. trial distance.
%
% Usage:
%   iic = IntertrialImageCorrelation(v, data_dir);
%   iic.neighborhood = 25;                     % optional, default 20
%
%   % optional score-based cell selection
%   scores = BaselineDriftAnalysis(v).getContributionScores('linear');
%   iic.scores = scores;
%   iic.score_quantile_range = [0, 0.2];       % keep lowest drift scores 20 % (high negative drift)
%
%   iic = iic.compute();
%   iic.plot();

    properties
        v                    ExperimentViewer
        data_dir             char
        neighborhood         double = 20     % pixel padding around ROI bounding box
        scores               cell   = {}     % {n_subj} of [N x 1] drift scores (optional)
        score_quantile_range double = [0 1]  % [q_lo, q_hi] on pooled score distribution
    end

    properties (SetAccess = private)
        corr_sheared_dff   double   % [2*T-1 x T x N_total]
        corr_sheared_anat  double   % [2*T-1 x T x N_total]
        n_trials           double   % T — set from first subject; others must match
    end

    % ------------------------------------------------------------------
    methods

        function obj = IntertrialImageCorrelation(v, data_dir, options)
            arguments
                v         ExperimentViewer
                data_dir  char
                options.neighborhood double = 20
            end
            obj.v            = v;
            obj.data_dir     = data_dir;
            obj.neighborhood = options.neighborhood;
        end

        % --------------------------------------------------------------
        function obj = compute(obj)
        % compute  Run pipeline steps 1–5; results in corr_sheared_dff/anat.

            traces_list = obj.v.filtered_traces;
            n_subj      = numel(traces_list);

            % Build per-subject cell inclusion masks from scores (once).
            cell_masks = obj.buildCellMasks(traces_list);

            all_sheared_dff  = [];
            all_sheared_anat = [];
            obj.n_trials     = [];

            for si = 1 : n_subj
                traces  = traces_list{si};
                subj_id = traces.subject_locations.subject_ID;
                n_incl  = sum(cell_masks{si});
                fprintf('[%d/%d] %s  (%d / %d cells)\n', ...
                    si, n_subj, subj_id, n_incl, traces.N);

                % -- (1a) Load dF/F stack ---------------------------------
                dff_files = dir(fullfile(obj.data_dir, [subj_id, '*.mat']));
                if isempty(dff_files)
                    warning('IntertrialImageCorrelation:dffNotFound', ...
                        'No dF/F file matching "%s*" in %s — skipping subject.', ...
                        subj_id, obj.data_dir);
                    continue
                end
                dff_stack = load( ...
                    fullfile(dff_files(1).folder, dff_files(1).name), 'stack').stack;
                [H, W, T] = size(dff_stack);

                if isempty(obj.n_trials)
                    obj.n_trials = T;
                elseif T ~= obj.n_trials
                    warning('IntertrialImageCorrelation:nTrialsMismatch', ...
                        '%s has %d trials (expected %d) — skipping subject.', ...
                        subj_id, T, obj.n_trials);
                    continue
                end

                % -- (1b) Load anatomy stack ------------------------------
                fish_path  = fullfile( ...
                    traces.subject_locations.subject_datapath, 'fish1.mat');
                anat_stack = loadAnatomyStack(fish_path, T);

                % -- (2–4) Crop ROIs, correlate, shear --------------------
                sheared_dff = shearStack( ...
                    IntertrialImageCorrelation.cropAndCorrelate( ...
                        dff_stack, traces, obj.neighborhood, cell_masks{si}), T);

                all_sheared_dff = cat(3, all_sheared_dff, sheared_dff);

                if ~isempty(anat_stack)
                    sheared_anat = shearStack( ...
                        IntertrialImageCorrelation.cropAndCorrelate( ...
                            anat_stack, traces, obj.neighborhood, cell_masks{si}), T);
                    all_sheared_anat = cat(3, all_sheared_anat, sheared_anat);
                end
            end

            obj.corr_sheared_dff  = all_sheared_dff;
            obj.corr_sheared_anat = all_sheared_anat;
            fprintf('Done. N_total = %d cells.\n', size(all_sheared_dff, 3));
        end

        % --------------------------------------------------------------
        function [hf, ax] = plot(obj, options)
        % plot  Two-curve line + shade plot: dF/F and anatomy (step 6).
            arguments
                obj
                options.color_dff   double    = [0.20 0.40 0.80]
                options.color_anat  double    = [0.60 0.60 0.60]
                options.cfg         PlotConfig = PlotConfig()
            end

            T         = obj.n_trials;
            distances = (-(T-1) : (T-1))';

            hf = figure('Color', 'w');
            ax = axes(hf);
            hold(ax, 'on');

            b(1) = plotCurve(ax, obj.corr_sheared_dff,  T, distances, ...
                options.color_dff,  options.cfg);
            b(2) = plotCurve(ax, obj.corr_sheared_anat, T, distances, ...
                options.color_anat, options.cfg);

            xline(ax, 0, '--k', 'Alpha', 0.4);
            xlabel(ax, 'Trial distance');
            ylabel(ax, 'Correlation (r)');
            legend(b, {'dF/F', 'Anatomy'}, 'Location', 'best');
            box(ax, 'off');
        end

    end % public methods

    % ------------------------------------------------------------------
    methods (Access = private)

        function masks = buildCellMasks(obj, traces_list)
        % buildCellMasks  Build per-subject logical inclusion masks.
        %
        %   If obj.scores is empty, all cells are included.
        %   Otherwise scores are pooled across subjects, quantile bounds
        %   are computed from obj.score_quantile_range, and each cell is
        %   included iff its score falls within [q_lo_val, q_hi_val].
        %
        %   Cells with non-finite scores are excluded when scores are
        %   provided, included otherwise.

            n_subj = numel(traces_list);
            masks  = cellfun(@(t) true(t.N, 1), traces_list, ...
                             'UniformOutput', false);

            if isempty(obj.scores); return; end

            % Pool all finite scores to compute global quantile bounds.
            all_scores = cell2mat( ...
                cellfun(@(s) s(:), obj.scores, 'UniformOutput', false));
            all_scores = all_scores(isfinite(all_scores));

            q_lo_val = quantile(all_scores, obj.score_quantile_range(1));
            q_hi_val = quantile(all_scores, obj.score_quantile_range(2));

            fprintf('Score filter: [%.3f, %.3f]  (quantile range [%.2f, %.2f])\n', ...
                q_lo_val, q_hi_val, ...
                obj.score_quantile_range(1), obj.score_quantile_range(2));

            for si = 1 : min(n_subj, numel(obj.scores))
                s = obj.scores{si}(:);
                N = traces_list{si}.N;
                if numel(s) ~= N
                    warning('IntertrialImageCorrelation:scoresSizeMismatch', ...
                        'scores{%d} has %d entries but traces has N=%d — ignoring scores for this subject.', ...
                        si, numel(s), N);
                    continue
                end
                masks{si} = isfinite(s) & s >= q_lo_val & s <= q_hi_val;
            end
        end

    end % private methods

    % ------------------------------------------------------------------
    methods (Static, Access = private)

        function roi_corr = cropAndCorrelate(stack, traces, neighborhood, cell_mask)
        % cropAndCorrelate  Steps 2–3: crop each ROI and compute [T x T x N].
        %
        %   stack        [H x W x T] image stack
        %   traces       ActivityTraces with ROImap and goodNeuron_IDs
        %   neighborhood pixel padding around ROI bounding box
        %   cell_mask    logical [N x 1]; excluded cells keep NaN slice
        %
        %   roi_corr     [T x T x N] intertrial Pearson correlation matrices

            [H, W, T] = size(stack);
            N         = traces.N;
            roi_corr  = nan(T, T, N);

            % Extract plain arrays before parfor — workers cannot broadcast
            % arbitrary class objects.
            goodNeuron_IDs = traces.goodNeuron_IDs;
            ROImap         = traces.ROImap;

            parfor ci = 1 : N
                if ~cell_mask(ci); continue; end

                roi_id        = goodNeuron_IDs(ci);
                roi_mask_full = (ROImap == roi_id);
                if ~any(roi_mask_full(:)); continue; end

                % Inline cropBounds: file-private functions are not reliably
                % accessible inside parfor workers.
                row_coords = find(any(roi_mask_full, 2));
                col_coords = find(any(roi_mask_full, 1));
                r1 = max(1, min(row_coords) - neighborhood);
                r2 = min(H, max(row_coords) + neighborhood);
                c1 = max(1, min(col_coords) - neighborhood);
                c2 = min(W, max(col_coords) + neighborhood);

                crop   = stack(r1:r2, c1:c2, :);              % [H_c x W_c x T]
                pixels = reshape(crop, (r2-r1+1)*(c2-c1+1), T);  % [n_px x T]
                roi_corr(:, :, ci) = corr(pixels, 'rows', 'complete');
            end
        end

    end % private static methods

end % classdef


% ======================================================================
% File-private helpers
% ======================================================================

function anat_stack = loadAnatomyStack(fish_path, expected_T)
% Load anatomy_imgs from the subject file; return [] if unavailable or
% if the trial count does not match expected_T.
    anat_stack = [];
    if ~exist(fish_path, 'file')
        warning('IntertrialImageCorrelation:subjectFileNotFound', ...
            'Subject file not found — anatomy skipped:\n  %s', fish_path);
        return
    end
    subject = robust_io('load', fish_path).fish1;
    anat    = subject.anatomy_imgs;
    if isempty(anat)
        warning('IntertrialImageCorrelation:anatomyEmpty', ...
            'anatomy_imgs is empty in %s — anatomy skipped.', fish_path);
        return
    end
    if size(anat, 3) ~= expected_T
        warning('IntertrialImageCorrelation:anatomyTrialsMismatch', ...
            'anatomy_imgs has %d trials (expected %d) — anatomy skipped.', ...
            size(anat, 3), expected_T);
        return
    end
    anat_stack = double(anat);
end

% ----------------------------------------------------------------------

function [r1, r2, c1, c2] = cropBounds(roi_mask_full, H, W, nb)
% Bounding box of roi_mask_full expanded by nb pixels, clamped to [H x W].
    row_coords = find(any(roi_mask_full, 2));
    col_coords = find(any(roi_mask_full, 1));
    r1 = max(1, min(row_coords) - nb);
    r2 = min(H, max(row_coords) + nb);
    c1 = max(1, min(col_coords) - nb);
    c2 = min(W, max(col_coords) + nb);
end

% ----------------------------------------------------------------------

function out = shearStack(C, T)
% Shear [T x T x N] into [2T-1 x T x N]: column i shifted up by (i-1)
% so the diagonal lands at row T.
    N   = size(C, 3);
    out = nan(2*T-1, T, N);
    for i = 1 : T
        row_start = T + 1 - i;
        out(row_start : row_start+T-1, i, :) = C(:, i, :);
    end
end

% ----------------------------------------------------------------------

function h = plotCurve(ax, S, T, distances, color, cfg)
% Compute mean ± SEM over the 3rd dim of sheared stack S and plot.
    if isempty(S); return; end
    n_dist = 2*T - 1;
    mu     = nan(n_dist, 1);
    sem    = nan(n_dist, 1);
    for r = 1 : n_dist
        vals = squeeze(S(r, :, :));
        vals = vals(isfinite(vals(:)));
        if numel(vals) < 2; continue; end
        mu(r)  = mean(vals);
        sem(r) = std(vals) / sqrt(numel(vals));
    end
    h = plotLineNShade(distances', mu', sem', color, cfg);
end
