classdef Behavior2PPlotter
% BEHAVIOR2PPLOTTER  Visualization and analysis tools for behavior-2p data.
%
%   Provides tail-motion onset/offset detection, peri-event neural PSTHs,
%   and brain-activity vs. breathing-rate regression for experiments where
%   calcium imaging and behavioural tracking are co-registered at 2-p
%   framerate.
%
%   PROPERTIES
%     v        - ExperimentViewer used to retrieve all data
%     n_std    - threshold multiplier: thr = mean + n_std * std  (default 2)
%     buffer_s - gap-fill debouncing window in seconds          (default 1)
%
%   EXAMPLE USAGE
%     bpp = Behavior2PPlotter(v);
%     bpp.n_std    = 2;
%     bpp.buffer_s = 1;
%
%     v.dataFilter.traceType = 'dFoverF_good';
%     v.dataFilter.interval  = [-5 50];
%
%     [onsets, thr] = bpp.getTailMotionOnsets();
%     hf            = bpp.plotPSTH(onsets, [-2 5]);
%     hf            = bpp.regressByBreathingRate();
%     hf            = bpp.plotBreathingEventPSTH([-0.5 1.5]);
%     hf            = bpp.plotBreathingCyclePSTH(100);
%
%   NOTE ON THRESHOLD
%     Thresholds are always computed from the FULL tail_motion_2p trace
%     (entire session, all trials) so that the same threshold is applied
%     regardless of which interval window is currently set in dataFilter.
%
%   SEE ALSO: ExperimentViewer, DataFilter, ModeSelector

    properties
        v        ExperimentViewer   % source data viewer
        n_std    double = 2         % threshold = mean + n_std * std
        buffer_s double = 1         % gap-fill debouncing buffer (seconds)
    end

    % =======================================================================
    methods

        function obj = Behavior2PPlotter(v)
            arguments
                v ExperimentViewer
            end
            obj.v = v;
        end

        % -------------------------------------------------------------------
        %  TAIL-MOTION EVENT DETECTION
        % -------------------------------------------------------------------

        function [onsets, thresholds] = getTailMotionOnsets(obj)
        % GETTAILMOTIONONSETS  Detect tail-motion onset times per fish.
        %
        %   Rising-edge crossings of the per-fish threshold are detected on
        %   the tail_motion_2p trace at the current dataFilter.interval.
        %   A gap-fill buffer of buffer_s seconds prevents brief dips below
        %   threshold from splitting a single bout into multiple events.
        %
        %   OUTPUTS
        %     onsets     {n_subjects x 1} cell of structs with fields:
        %                  .trial_idx  [n_events x 1]  trial number
        %                  .time_s     [n_events x 1]  onset time (s) within
        %                              the current peri-event window
        %     thresholds [n_subjects x 1] threshold value per fish
            [onsets, thresholds] = obj.detectTailMotionEvents('onset');
        end

        function [offsets, thresholds] = getTailMotionOffsets(obj)
        % GETTAILMOTIONOFFSETS  Detect tail-motion offset times per fish.
        %
        %   Symmetric to getTailMotionOnsets; detects falling-edge crossings.
        %   The gap-fill buffer prevents brief spikes above threshold from
        %   splitting a single rest period into multiple events.
        %
        %   OUTPUTS: see getTailMotionOnsets
            [offsets, thresholds] = obj.detectTailMotionEvents('offset');
        end

        % -------------------------------------------------------------------
        %  BEHAVIOR IMAGESC
        % -------------------------------------------------------------------

        function hf = plotBehaviorImagesc(obj, varargin)
        % PLOTBEHAVIORIMAGESC  Imagesc of a behavior2p trace across trials.
        %
        %   Plots the trace specified by v.dataFilter.traceType as a 2-D
        %   imagesc (time × trials), averaged across subjects.  Mirrors the
        %   printImagesc / prepTraces pattern from the analysis scripts.
        %
        %   Set v.dataFilter.traceType before calling:
        %     'tail_motion_2p'         tail motion at 2p framerate
        %     'breathing_inst_freq_2p' instantaneous breathing frequency
        %     (or any other behavior2p tag)
        %
        %   NAME-VALUE ARGUMENTS
        %     'dozscore'  logical  z-score each subject using its own
        %                          mean and std computed on the FULL trace
        %                          (interval = []).  Default: true.
        %
        %   OUTPUT
        %     hf  figure handle

            % --- parse varargin -------------------------------------------
            dozscore = true;
            for k = 1:2:numel(varargin)
                if strcmpi(varargin{k}, 'dozscore')
                    dozscore = varargin{k+1};
                end
            end

            v      = obj.v;
            cfg    = v.plotConfig;
            ps_lim = v.dataFilter.interval;

            % --- per-subject z-score statistics from full trace -----------
            if dozscore
                dft_full = v.dataFilter;
                dft_full.interval = [];
                v_full = v;
                v_full.dataFilter = dft_full;
                [ev_full, ~] = dft_full.filterData(v_full);
                has_full = ~cellfun(@isempty, ev_full);
                ev_full  = ev_full(has_full);
                mu  = cellfun(@(x) mean(x(:), 'omitmissing'), ev_full);
                err = cellfun(@(x) std(x(:),  'omitmissing'), ev_full);
            end

            % --- traces at current interval --------------------------------
            [events, labs] = v.dataFilter.filterData(v);
            hasdata = ~cellfun(@isempty, events);
            events  = events(hasdata);
            labs    = labs(hasdata);
            events  = cellfun(@squeeze, events, 'UniformOutput', false);

            if isempty(events)
                warning('Behavior2PPlotter:plotBehaviorImagesc: no data.');
                hf = []; return;
            end

            % --- identity normalization when z-score is off ---------------
            if ~dozscore
                nsubj = numel(events);
                mu    = zeros(nsubj, 1);
                err   = ones(nsubj, 1);
            end

            [L, ntrials] = size(events{1});
            nsubj = numel(events);
            if isempty(ps_lim)
                t = (0:L-1) / v.filtered_traces{1}.framerate;
            else
                t = linspace(ps_lim(1), ps_lim(2), L);
            end

            % --- build [L x ntrials x nsubj] and z-score -----------------
            M = nan(L, ntrials, nsubj);
            for i = 1:nsubj
                M(:,:,i) = (events{i} - mu(i)) ./ err(i);
            end

            avg = mean(M, 3, 'omitmissing');   % [L x ntrials]

            % --- plot ------------------------------------------------------
            hf = figure;
            imagesc(1:ntrials, t, avg);
            xticks(1:ntrials);
            xticklabels(labs{1});
            axis square;
            colormap(cfg.colormapName);
            xlabel('Trial');
            ylabel('Time (s)');
            cb = colorbar('Color', cfg.axcol);
            if dozscore
                cb.Label.String = [v.dataFilter.traceType, ' (z-score)'];
            else
                cb.Label.String = v.dataFilter.traceType;
            end
            cb.Label.FontSize = gca().FontSize;
            set(gca, 'color', cfg.bgcol, ...
                'XColor', cfg.axcol, 'YColor', cfg.axcol, 'ZColor', cfg.axcol);
            set(gcf, 'color', cfg.bgcol);
        end

        % -------------------------------------------------------------------
        %  NEURAL ACTIVITY LOADER
        % -------------------------------------------------------------------

        function [neural_data, t] = loadNeuralActivity(obj)
        % LOADNEURALACTIVITY  Load peri-event neural activity averaged over units.
        %
        %   Uses ModeSelector with the current dataFilter settings (traceType,
        %   interval, mode_name, etc.) to extract [T x N x n_trials] data,
        %   then averages over the N (units/modes) dimension.
        %
        %   Tip: set v.dataFilter.traceType to your preferred neural trace
        %   (e.g., 'dFoverF_good', 'pSpike') before calling.
        %
        %   OUTPUTS
        %     neural_data {n_subjects x 1} cell of [T x n_trials] double
        %     t           [1 x T] time axis in seconds

            v     = obj.v;
            ps_lim = v.dataFilter.interval;

            [~, events, ~] = ModeSelector(v).extract;
            n = numel(events);
            neural_data = cell(n, 1);

            for i = 1:n
                if isempty(events{i}); continue; end
                % events{i}: [T, N, n_trials]
                nd = mean(events{i}, 2, 'omitmissing');  % [T, 1, n_trials]
                nd = squeeze(nd);                          % [T, n_trials]
                if isvector(nd); nd = nd(:); end           % single trial
                neural_data{i} = nd;
            end

            % Time axis from first non-empty subject
            first = find(~cellfun(@isempty, neural_data), 1);
            if isempty(first); t = []; return; end
            T = size(neural_data{first}, 1);
            if isempty(ps_lim)
                t = (0:T-1) / v.filtered_traces{first}.framerate;
            else
                t = linspace(ps_lim(1), ps_lim(2), T);
            end
        end

        % -------------------------------------------------------------------
        %  PERI-EVENT NEURAL PSTH
        % -------------------------------------------------------------------

        function [hf, snippets, t_psth] = plotPSTH(obj, events, trange)
        % PLOTPSTH  Plot peri-event neural PSTH with tail-motion underlay.
        %
        %   Cuts windows of mean neural activity (averaged over units) and the
        %   corresponding tail_motion_2p signal around each detected event.
        %   Mean ± SEM across all pooled events is plotted for both signals
        %   on a dual y-axis so that the independent scales are preserved.
        %
        %   OPTIONAL OUTPUTS
        %     snippets  [n_psth_frames x n_events] raw neural snippets
        %               (one column per event, pooled across all subjects)
        %               Suitable for passing to extractActivityMetric as
        %               'motion tuning' with 'PreEventFrames'.
        %     t_psth    [1 x n_psth_frames] time axis in seconds
        %
        %   INPUTS
        %     events  cell output from getTailMotionOnsets / getTailMotionOffsets
        %     trange  [pre_s, post_s] window around each event (default [-2 5])
        %
        %   OUTPUT
        %     hf figure handle
            arguments
                obj
                events cell
                trange (1,2) double = [-2, 5]
            end

            v      = obj.v;
            cfg    = v.plotConfig;
            ps_lim = v.dataFilter.interval;
            traces = v.filtered_traces;
            n      = numel(events);

            % Load neural activity [T x n_trials] per subject
            [neural_data, ~] = obj.loadNeuralActivity();

            % Load tail_motion_2p at the same peri-event interval
            dft_beh = v.dataFilter;
            dft_beh.traceType = 'tail_motion_2p';
            v_beh = v;
            v_beh.dataFilter = dft_beh;
            [beh_events, ~] = dft_beh.filterData(v_beh);
            % beh_events{i}: [T_win, 1, n_trials]

            tail_data = cell(n, 1);
            for i = 1:n
                if isempty(beh_events{i}); continue; end
                td = squeeze(beh_events{i});     % [T_win, n_trials]
                if isvector(td); td = td(:); end
                tail_data{i} = td;
            end

            % Reference framerate from first usable subject
            fs_ref = [];
            for i = 1:n
                if ~isempty(neural_data{i})
                    fs_ref = traces{i}.framerate;
                    break;
                end
            end
            if isempty(fs_ref)
                warning('Behavior2PPlotter:plotPSTH: no neural data available.');
                hf = [];
                return;
            end

            n_psth_fr     = floor(diff(trange) * fs_ref);
            t_psth        = trange(1) + (0:n_psth_fr-1) / fs_ref;
            snippets      = nan(n_psth_fr, 0);   % neural
            tail_snippets = nan(n_psth_fr, 0);   % tail motion

            for i = 1:n
                if isempty(events{i}) || isempty(neural_data{i}); continue; end

                ev     = events{i};
                nd     = neural_data{i};    % [T_neural, n_trials]
                td     = tail_data{i};      % [T_tail,   n_trials] or []
                fs     = traces{i}.framerate;
                n_fr_i = floor(diff(trange) * fs);

                for k = 1:numel(ev.time_s)
                    t_on = ev.time_s(k);
                    j    = ev.trial_idx(k);
                    if j > size(nd, 2); continue; end

                    % Convert event time to frame index
                    if isempty(ps_lim)
                        f_on = round(t_on * fs) + 1;
                    else
                        f_on = round((t_on - ps_lim(1)) * fs) + 1;
                    end
                    f_start = f_on + round(trange(1) * fs);
                    f_end   = f_start + n_fr_i - 1;

                    if f_start < 1 || f_end > size(nd, 1); continue; end

                    % Neural snippet
                    snippet = nd(f_start:f_end, j);
                    if numel(snippet) ~= n_psth_fr
                        snippet = interp1(linspace(0,1,numel(snippet)), ...
                            double(snippet), linspace(0,1,n_psth_fr))';
                    end
                    snippets(:, end+1) = snippet(:); %#ok<AGROW>

                    % Tail-motion snippet (same window; NaN-pad if unavailable)
                    if ~isempty(td) && j <= size(td,2) && f_end <= size(td,1)
                        tail_snip = td(f_start:f_end, j);
                        if numel(tail_snip) ~= n_psth_fr
                            tail_snip = interp1(linspace(0,1,numel(tail_snip)), ...
                                double(tail_snip), linspace(0,1,n_psth_fr))';
                        end
                    else
                        tail_snip = nan(n_psth_fr, 1);
                    end
                    tail_snippets(:, end+1) = tail_snip(:); %#ok<AGROW>
                end
            end

            if isempty(snippets)
                warning('Behavior2PPlotter:plotPSTH: no events within bounds.');
                hf = [];
                return;
            end

            % Compute mean ± SEM for each signal
            n_ev  = size(snippets, 2);
            mu_n  = mean(snippets,      2, 'omitmissing');
            err_n = std(snippets,  0,   2, 'omitmissing') / sqrt(n_ev);
            mu_t  = mean(tail_snippets, 2, 'omitmissing');
            n_t   = sum(~isnan(tail_snippets), 2);   % valid count per timepoint
            err_t = std(tail_snippets, 0, 2, 'omitmissing') ./ sqrt(max(n_t, 1));

            col_neural = cfg.c(1,:);
            col_tail   = cfg.c(2,:);

            hf = figure;

            % Left y-axis: neural activity
            yyaxis left;
            plotLineNShade(t_psth, mu_n, err_n, col_neural, cfg);
            hold on;
            ylabel('Mean neural activity (a.u.)');
            set(gca, 'YColor', col_neural);

            % Right y-axis: tail motion
            yyaxis right;
            hold on;
            plotLineNShade(t_psth, mu_t, err_t, col_tail, cfg);
            ylabel('Tail motion (a.u.)');
            set(gca, 'YColor', col_tail);

            % Shared decorations (drawn on whichever side is active — cosmetic only)
            xline(0, '--', 'Color', cfg.axcol, 'LineWidth', 1);
            xlabel('Time from tail-motion onset (s)');
            title(sprintf('n = %d events', n_ev), 'Color', cfg.textcol);
            axis tight; box off;
            set(gca, 'color', cfg.bgcol, 'XColor', cfg.axcol);
            set(gcf, 'color', cfg.bgcol);
        end

        % -------------------------------------------------------------------
        %  PER-UNIT MOTION SNIPPETS
        % -------------------------------------------------------------------

        function [snippets, t_psth, pre_event_frames] = getMotionSnippets(obj, events, trange)
        % GETMOTIONSNIPPETS  Extract per-unit peri-event neural snippets.
        %
        %   Unlike plotPSTH (which collapses over units for visualisation),
        %   this method preserves the unit dimension and returns one snippet
        %   per unit per event, suitable for per-unit 'motion tuning' scoring
        %   via extractActivityMetric.
        %
        %   INPUTS
        %     events  cell from getTailMotionOnsets / getTailMotionOffsets
        %     trange  [pre_s, post_s] window in seconds (default [-2 5])
        %
        %   OUTPUTS
        %     snippets         {n_subjects} cell of [n_frames x n_units x n_events]
        %     t_psth           [1 x n_frames] common time axis (seconds)
        %     pre_event_frames integer — number of frames before t = 0;
        %                      pass directly as 'PreEventFrames' to
        %                      extractActivityMetric('motion tuning', ...).
        %
        %   EXAMPLE
        %     [onsets, ~]        = bpp.getTailMotionOnsets();
        %     [offsets, ~]       = bpp.getTailMotionOffsets();
        %     [snips, ~, pre_fr] = bpp.getMotionSnippets(onsets, [-2 5]);
        %     [~, full_nd, ~]    = ModeSelector(bpp.v).extract;
        %     ps = bpp.v.dataFilter.interval;
        %     fs = bpp.v.filtered_traces{i}.framerate;
        %     for i = 1:numel(snips)
        %         tuning{i} = extractActivityMetric(snips{i}, ...
        %             'motion tuning', 'cells', ...
        %             'PreEventFrames', pre_fr, ...
        %             'FullData',       full_nd{i}, ...
        %             'EventOnsets',    onsets{i}, ...
        %             'EventOffsets',   offsets{i}, ...
        %             'FrameRate',      fs, ...
        %             'PsLim',          ps);
        %     end
            arguments
                obj
                events cell
                trange (1,2) double = [-2, 5]
            end

            v      = obj.v;
            ps_lim = v.dataFilter.interval;
            traces = v.filtered_traces;
            n      = numel(events);

            % Raw neural data [T x N x n_trials] per subject (no unit averaging)
            [~, raw_data, ~] = ModeSelector(v).extract;

            % Reference framerate from first usable subject
            fs_ref = [];
            for i = 1:n
                if ~isempty(raw_data{i})
                    fs_ref = traces{i}.framerate;
                    break;
                end
            end
            if isempty(fs_ref)
                warning('Behavior2PPlotter:getMotionSnippets: no neural data.');
                snippets = {}; t_psth = []; pre_event_frames = 0;
                return;
            end

            n_psth_fr        = floor(diff(trange) * fs_ref);
            t_psth           = trange(1) + (0:n_psth_fr-1) / fs_ref;
            pre_event_frames = floor(abs(trange(1)) * fs_ref);

            snippets = cell(n, 1);

            for i = 1:n
                if isempty(events{i}) || isempty(raw_data{i}); continue; end

                ev     = events{i};
                nd     = raw_data{i};    % [T x N x n_trials]
                fs     = traces{i}.framerate;
                N      = size(nd, 2);
                n_fr_i = floor(diff(trange) * fs);

                snips_i = nan(n_psth_fr, N, 0);

                for k = 1:numel(ev.time_s)
                    t_on = ev.time_s(k);
                    j    = ev.trial_idx(k);
                    if j > size(nd, 3); continue; end

                    trial_nd = nd(:, :, j);   % [T x N]

                    if isempty(ps_lim)
                        f_on = round(t_on * fs) + 1;
                    else
                        f_on = round((t_on - ps_lim(1)) * fs) + 1;
                    end
                    f_start = f_on + round(trange(1) * fs);
                    f_end   = f_start + n_fr_i - 1;

                    if f_start < 1 || f_end > size(nd, 1); continue; end

                    snippet = trial_nd(f_start:f_end, :);   % [n_fr_i x N]
                    if n_fr_i ~= n_psth_fr
                        % interp1 on matrix Y interpolates each column
                        snippet = interp1(linspace(0,1,n_fr_i), ...
                            double(snippet), linspace(0,1,n_psth_fr));
                    end
                    snips_i(:, :, end+1) = snippet;   %#ok<AGROW>
                end

                snippets{i} = snips_i;   % [n_psth_fr x N x n_events_i]
            end
        end

        % -------------------------------------------------------------------
        %  BRAIN ACTIVITY vs. BREATHING RATE REGRESSION
        % -------------------------------------------------------------------

        function [hf, lme] = regressByBreathingRate(obj)
        % REGRESSBYBREATHINGRATE  Regress instantaneous neural activity on breathing rate.
        %
        %   For every time point within every trial, pairs the mean neural
        %   activity (averaged over units) with the simultaneously recorded
        %   instantaneous breathing rate (breathing_inst_freq_2p).  Both
        %   signals are at the 2-p framerate and are aligned via DataFilter
        %   to the same peri-event window.
        %
        %   Fits a linear mixed-effects model accounting for the nested
        %   structure (time points within subjects):
        %
        %       brain_activity ~ 1 + breathing_rate + (1|subject)
        %
        %   OUTPUTS
        %     hf  figure handle
        %     lme fitted LinearMixedModel object (from fitlme)
        %
        %   PLOT
        %     Scatter of all (breathing_rate, brain_activity) pairs coloured
        %     by subject, with the LME fixed-effect trend line overlaid.
        %     Title reports slope, significance stars, and marginal R².

            v      = obj.v;
            cfg    = v.plotConfig;
            traces = v.filtered_traces;
            n      = numel(traces);

            % --- Load neural activity [T x n_trials] per subject ----------
            [neural_data, ~] = obj.loadNeuralActivity();

            % --- Load breathing rate at the same peri-event window --------
            dft = v.dataFilter;
            dft.traceType = 'breathing_inst_freq_2p';
            v_beh = v;
            v_beh.dataFilter = dft;
            [beh_events, ~] = dft.filterData(v_beh);
            % beh_events{i}: [T, 1, n_trials]

            % --- Accumulate pooled data for LME ---------------------------
            all_y    = [];   % mean neural activity [N_pts x 1]
            all_x    = [];   % breathing rate       [N_pts x 1]
            all_subj = {};   % subject label        [N_pts x 1]

            subj_labels = cell(n, 1);
            subj_colors = cell(n, 1);

            for i = 1:n
                nd  = neural_data{i};    % [T, n_trials]
                beh = beh_events{i};     % [T, 1, n_trials]
                if isempty(nd) || isempty(beh); continue; end

                beh = squeeze(beh);      % [T, n_trials]
                if isvector(beh); beh = beh(:); end

                T_use    = min(size(nd, 1), size(beh, 1));
                ntr_use  = min(size(nd, 2), size(beh, 2));
                subj_id  = traces{i}.subject_locations.subject_ID;
                subj_labels{i} = subj_id;
                subj_colors{i} = cfg.c(mod(i-1, size(cfg.c,1))+1, :);

                for j = 1:ntr_use
                    y_tr  = nd(1:T_use, j);
                    x_tr  = beh(1:T_use, j);
                    valid = ~isnan(y_tr) & ~isnan(x_tr);
                    if ~any(valid); continue; end
                    nv = sum(valid);
                    all_y    = [all_y;    y_tr(valid)];   %#ok<AGROW>
                    all_x    = [all_x;    x_tr(valid)];   %#ok<AGROW>
                    all_subj = [all_subj; repmat({subj_id}, nv, 1)]; %#ok<AGROW>
                end
            end

            if isempty(all_y)
                warning('Behavior2PPlotter:regressByBreathingRate: no valid data.');
                hf = []; lme = [];
                return;
            end

            % --- Fit LME: brain_activity ~ 1 + breathing_rate + (1|subj) -
            lme_result = statsUtils.lmeRegress(all_y, all_x, all_subj);
            lme     = lme_result.lme;
            slope   = lme_result.slope;
            slope_p = lme_result.p;
            r2_marg = lme_result.r2_marginal;

            % --- Plot scatter + trend line --------------------------------
            hf = figure;
            hold on;
            valid_subj = ~cellfun(@isempty, subj_labels);
            leg_handles = gobjects(sum(valid_subj), 1);
            leg_idx = 0;
            for i = 1:n
                if isempty(subj_labels{i}); continue; end
                idx = strcmp(all_subj, subj_labels{i});
                leg_idx = leg_idx + 1;
                leg_handles(leg_idx) = scatter( ...
                    all_x(idx), all_y(idx), ...
                    8, subj_colors{i}, 'filled', ...
                    'MarkerFaceAlpha', 0.3, ...
                    'DisplayName', subj_labels{i});
            end

            % Fixed-effects trend line
            x_rng = linspace(min(all_x), max(all_x), 200);
            y_rng = lme_result.intercept + lme_result.slope * x_rng;
            plot(x_rng, y_rng, '-', 'Color', cfg.axcol, 'LineWidth', 2, ...
                'DisplayName', 'LME fit');

            xlabel('Breathing rate (Hz)');
            ylabel('Mean neural activity (a.u.)');
            sig_str = statsUtils.pvalToStars(slope_p);
            title(sprintf('slope = %.3g,  %s,  R²_{marg} = %.3f', ...
                slope, sig_str, r2_marg), 'Color', cfg.textcol);
            legend(leg_handles, 'Location', 'best', ...
                'TextColor', cfg.textcol, 'Color', cfg.bgcol);
            axis tight; box off;
            set(gca, 'color', cfg.bgcol, ...
                'XColor', cfg.axcol, 'YColor', cfg.axcol, 'ZColor', cfg.axcol);
            set(gcf, 'color', cfg.bgcol);
        end

    end % public methods

    % =======================================================================
    methods (Access = private)

        function thresholds = getThresholds(obj)
        % GETTHRESHOLDS  Per-subject tail-motion threshold from full session.
        %
        %   Computed on the FULL tail_motion_2p trace via getBehavior2PTrace,
        %   ignoring any interval currently set in dataFilter, so the same
        %   threshold applies to all interval windows.
            traces = obj.v.filtered_traces;
            n = numel(traces);
            thresholds = nan(n, 1);
            for i = 1:n
                [raw, ~] = traces{i}.getBehavior2PTrace('tail_motion_2p');
                if isempty(raw); continue; end
                x = raw(:);
                thresholds(i) = mean(x, 'omitmissing') + ...
                    obj.n_std * std(x, 'omitmissing');
            end
        end

        function [events, thresholds] = detectTailMotionEvents(obj, direction)
        % DETECTTAILMOTIONEVENTS  Shared onset/offset detection logic.
            v          = obj.v;
            thresholds = obj.getThresholds();
            ps_lim     = v.dataFilter.interval;

            % Get peri-event tail_motion_2p at the current interval
            dft = v.dataFilter;
            dft.traceType = 'tail_motion_2p';
            v_beh = v;
            v_beh.dataFilter = dft;
            [beh_events, ~] = dft.filterData(v_beh);
            % beh_events{i}: [T_win, 1, n_trials]

            n      = numel(beh_events);
            events = cell(n, 1);

            for i = 1:n
                if isempty(beh_events{i}) || isnan(thresholds(i)); continue; end

                M = squeeze(beh_events{i});      % [T_win, n_trials]
                if isvector(M); M = M(:); end    % single-trial safety
                [T_win, n_trials] = size(M);

                fs            = v.filtered_traces{i}.framerate;
                buffer_frames = max(1, round(obj.buffer_s * fs));

                if isempty(ps_lim)
                    t_win = (0:T_win-1) / fs;
                else
                    t_win = linspace(ps_lim(1), ps_lim(2), T_win);
                end

                all_trial = [];
                all_time  = [];
                for j = 1:n_trials
                    e_frames = detectCrossings( ...
                        M(:, j), thresholds(i), buffer_frames, direction);
                    if isempty(e_frames); continue; end
                    all_trial = [all_trial; repmat(j, numel(e_frames), 1)]; %#ok<AGROW>
                    all_time  = [all_time;  t_win(e_frames)'];              %#ok<AGROW>
                end
                events{i} = struct('trial_idx', all_trial, 'time_s', all_time);
            end
        end

    end % private methods
end % classdef

% ==========================================================================
%  FILE-PRIVATE HELPER FUNCTIONS
% ==========================================================================

function frames = detectCrossings(sig, threshold, buffer_frames, direction)
% DETECTCROSSINGS  Debounced threshold-crossing detector.
%
%   Applies a gap-fill buffer to suppress transient excursions before
%   detecting crossings:
%     'onset'  – fills 0-runs shorter than buffer_frames, then finds rising edges
%     'offset' – fills 1-runs shorter than buffer_frames, then finds falling edges
%
%   INPUTS
%     sig           [T x 1] signal
%     threshold     scalar threshold
%     buffer_frames integer minimum run length to preserve
%     direction     'onset' or 'offset'
%
%   OUTPUT
%     frames integer vector of crossing frame indices

    above = double(sig(:) >= threshold);

    switch direction
        case 'onset'
            above  = fillSmallGaps(above, buffer_frames, 0);  % fill 0-gaps
            frames = find(diff([0; above]) == 1);
        case 'offset'
            above  = fillSmallGaps(above, buffer_frames, 1);  % fill 1-gaps
            frames = find(diff([above; 0]) == -1);
        otherwise
            error('Behavior2PPlotter:detectCrossings: unknown direction "%s"', direction);
    end
end

function sig = fillSmallGaps(sig, min_len, gap_val)
% FILLSMALLGAPS  Replace short runs of gap_val with (1 - gap_val).
%
%   Runs of gap_val strictly shorter than min_len are filled.  Longer runs
%   are left intact.
%
%   Example (gap_val=0, min_len=3):
%     input  [1 1 0 0 1 1 0 0 0 1]
%     output [1 1 1 1 1 1 0 0 0 1]   (2-frame gap filled, 3-frame gap kept)

    fill_val = 1 - gap_val;
    k = 1;
    n = numel(sig);
    while k <= n
        if sig(k) == gap_val
            j = k;
            while j <= n && sig(j) == gap_val
                j = j + 1;
            end
            if (j - k) < min_len
                sig(k:j-1) = fill_val;
            end
            k = j;
        else
            k = k + 1;
        end
    end
end
