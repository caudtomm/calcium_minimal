classdef SliceTCA_Analysis
% SliceTCA_Analysis  Slice TCA decomposition workflow for neural activity.
%
%   Implements the MATLAB side of a three-step sliceTCA pipeline:
%     1. Export per-subject data tensors to .mat files  (saveInputFiles)
%     2. [External] Run Python sliceTCA script on each file
%     3. Load, inspect and visualize results             (loadResults, plotSubject)
%
%   The data tensor written per subject follows sliceTCA convention:
%     [trials x neurons x time_bins]  (double)
%
%   Filtering (groups, stimuli, repetitions, time window) is read from
%   obj.v.dataFilter at call time.
%
%   Usage:
%       sla = SliceTCA_Analysis(v);
%       sla.output_dir = 'path/to/dir';
%
%       paths   = sla.saveInputFiles();          % step 1
%       % ... run Python script externally ...   % step 2
%       results = sla.loadResults();             % step 3a
%       hf      = sla.plotSubject(results, 1);   % step 3b

    properties
        v          ExperimentViewer
        output_dir char = ''   % default save location for saveInputFiles
    end

    methods

        function obj = SliceTCA_Analysis(v, output_dir)
            arguments
                v          ExperimentViewer
                output_dir char = ''
            end
            obj.v          = v;
            obj.output_dir = output_dir;
        end

        % ------------------------------------------------------------------
        function file_paths = saveInputFiles(obj, output_dir)
        % saveInputFiles  Extract filtered data and write per-subject .mat files.
        %
        %   file_paths = saveInputFiles()
        %   file_paths = saveInputFiles(output_dir)
        %
        %   Calls ModeSelector with the current v.dataFilter settings
        %   (mode, stimuli, interval, repetitions, etc.), forces
        %   chronological trial ordering, then writes one HDF5-compatible
        %   .mat file per subject with the following variables:
        %
        %     data          [trials x neurons x time]  double
        %     trial_labels  {trials x 1}  cell of char — stimulus identity
        %     t             [1 x time]    double       — time axis (seconds)
        %     subject_id    char                       — subject identifier
        %     group         char                       — subject group tag
        %     framerate     double                     — imaging rate (Hz)

            arguments
                obj        SliceTCA_Analysis
                output_dir char = obj.output_dir
            end

            if isempty(output_dir)
                error('SliceTCA_Analysis:noOutputDir', ...
                    ['Specify an output directory as argument or set ', ...
                     'obj.output_dir before calling saveInputFiles.']);
            end
            if ~isfolder(output_dir)
                mkdir(output_dir);
            end

            %% -- Extract filtered neural data ----------------------------
            % Force chronological ordering so trial index == recording order.
            prev_sorting = obj.v.dataFilter.trial_sorting;
            obj.v.dataFilter.trial_sorting = 'chronological';
            [~, events, labels] = ModeSelector(obj.v).extract;
            obj.v.dataFilter.trial_sorting = prev_sorting;

            n_subj = numel(events);
            ps_lim = obj.v.dataFilter.interval;   % [t_start, t_end] seconds
            if isempty(ps_lim)
                fs = obj.v.filtered_traces{1}.framerate;
                ps_lim = [0 size(events{1},1)/fs];
            end

            % Resolve subject IDs in the same order ModeSelector used
            [subj_ids, ~] = obj.v.dataFilter.getSubjectIDs(obj.v.subjectTab);
            file_paths    = cell(n_subj, 1);

            %% -- Write one file per subject ------------------------------
            for i = 1:n_subj
                data_tnt = events{i};     % [time x neurons x trials]
                if isempty(data_tnt)
                    warning('SliceTCA_Analysis:emptyData', ...
                        'Subject %d (%s): no data after filtering — skipping.', ...
                        i, subj_ids{i});
                    continue
                end

                % Permute to sliceTCA convention: [trials x neurons x time]
                data = permute(data_tnt, [3 2 1]);

                % Metadata
                trial_labels = labels{i};           % {trials x 1} stimulus strings
                n_time       = size(data, 3);
                traces       = obj.v.filtered_traces{i};
                framerate    = traces.framerate;
                t            = ps_lim(1) + (0:n_time - 1) / framerate;  % [1 x T] seconds
                subject_id   = subj_ids{i};
                group        = char(traces.subject_group);

                % Write as HDF5-compatible .mat (readable by Python mat73/h5py)
                safe_name  = matlab.lang.makeValidName(subject_id);
                fpath      = fullfiletol(output_dir, [safe_name, '_input.mat']);
                save(fpath, 'data', 'trial_labels', 't', ...
                     'subject_id', 'group', 'framerate', '-v7.3');

                file_paths{i} = fpath;
                fprintf('[SliceTCA] Saved %s  [%d trials x %d neurons x %d time bins]\n', ...
                    fpath, size(data,1), size(data,2), size(data,3));
            end
        end

        % ------------------------------------------------------------------
        function results = loadResults(obj, results_dir)
        % loadResults  Batch-load *_slicetca.mat result files.
        %
        %   results = loadResults(results_dir)
        %   results = loadResults()              % uses obj.output_dir
        %
        %   Returns a {n_files x 1} cell array of structs, one per subject,
        %   each holding the variables written by run_slicetca.py:
        %     subject_id, group, best_ranks, loss_grid, losses,
        %     reconstruction, t, framerate, trial_labels,
        %     components_{k}_scores   [r_k x sliced_dim]
        %     components_{k}_weights  [r_k x dim_j x dim_l]   k = 0,1,2

            arguments
                obj         SliceTCA_Analysis
                results_dir char = obj.output_dir
            end

            if isempty(results_dir)
                error('SliceTCA_Analysis:noResultsDir', ...
                    'Specify a results directory or set obj.output_dir.');
            end

            files = dir(fullfiletol(results_dir, '*_slicetca.mat'));
            if isempty(files)
                error('SliceTCA_Analysis:noResults', ...
                    'No *_slicetca.mat files found in: %s', results_dir);
            end

            n       = numel(files);
            results = cell(n, 1);
            for i = 1:n
                fpath      = fullfiletol(results_dir, files(i).name);
                results{i} = load(fpath);
                fprintf('[SliceTCA] Loaded [%d/%d]  %s\n', i, n, files(i).name);
            end
        end

        % ------------------------------------------------------------------
        function hf = plotSubject(obj, results, fish_idx) %#ok<INUSL>
        % plotSubject  Plot scores and weights for one subject.
        %
        %   hf = plotSubject(results, fish_idx)
        %
        %   results    {n_subj} cell from loadResults
        %   fish_idx   1-based index into results
        %
        %   Figure layout — 3 rows (partition) x 2 columns:
        %     Col 1: score matrix    [r_k x sliced_dim]   (imagesc)
        %     Col 2: weight heatmaps [r_k components stacked vertically]
        %
        %   Partition 0  trial-slice : scores (r0 x trials),   weights (r0 x neurons x time)
        %   Partition 1  neuron-slice: scores (r1 x neurons),  weights (r1 x trials  x time)
        %   Partition 2  time-slice  : scores (r2 x time),     weights (r2 x trials  x neurons)

            arguments
                obj      SliceTCA_Analysis
                results  cell
                fish_idx (1,1) double {mustBePositive, mustBeInteger}
            end

            if fish_idx > numel(results) || isempty(results{fish_idx})
                error('SliceTCA_Analysis:invalidIndex', ...
                    'fish_idx=%d out of range (results has %d entries).', ...
                    fish_idx, numel(results));
            end

            r   = results{fish_idx};
            t   = double(r.t(:)');          % [1 x S] time axis in seconds
            sid = '';
            if isfield(r, 'subject_id');  sid = char(r.subject_id);  end

            part_names   = {'Trial-slice',  'Neuron-slice', 'Time-slice'};
            score_xlbls  = {'Trial',        'Neuron',       'Time (s)'};
            weight_xlbls = {'Time (s)',     'Time (s)',     'Neuron'};
            weight_ylbls = {'Neuron',       'Trial',        'Trial'};
            % partition 2 scores use the time axis; others use integer index
            score_use_t  = [false, false, true];
            % partition 0 & 1 weight x-axes are time; partition 2 is neuron index
            weight_use_t = [true, true, false];

            hf = figure('Name', sprintf('SliceTCA — %s', sid), ...
                        'Color', 'w', 'Units', 'normalized', ...
                        'OuterPosition', [0 0 1 1]);
            tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

            for k = 0:2
                sk = sprintf('components_%d_scores',  k);
                wk = sprintf('components_%d_weights', k);
                if ~isfield(r, sk) || ~isfield(r, wk);  nexttile; nexttile; continue;  end

                scores  = double(r.(sk));    % [r_k x sliced_dim]
                weights = double(r.(wk));    % [r_k x dim_j x dim_l]

                r_k   = size(scores,  1);
                dim_j = size(weights, 2);
                dim_l = size(weights, 3);

                if r_k == 0;  nexttile; nexttile; continue;  end

                % ---- Scores ----------------------------------------------
                ax_s = nexttile;
                if score_use_t(k+1) && numel(t) == size(scores, 2)
                    imagesc(ax_s, t, 1:r_k, scores);
                else
                    imagesc(ax_s, scores);
                end
                colorbar(ax_s);
                xlabel(ax_s, score_xlbls{k+1});
                ylabel(ax_s, 'Component');
                yticks(ax_s, 1:r_k);
                title(ax_s, [part_names{k+1}, ' — scores']);
                axis(ax_s, 'tight');  box(ax_s, 'off');

                % ---- Weights ---------------------------------------------
                ax_w = nexttile;

                % Stack r_k weight matrices vertically: [r_k*dim_j x dim_l]
                % permute [r_k,dim_j,dim_l] → [dim_j,r_k,dim_l] then reshape
                stacked = reshape(permute(weights, [2 1 3]), dim_j * r_k, dim_l);

                if weight_use_t(k+1) && numel(t) == dim_l
                    imagesc(ax_w, t, 1:dim_j*r_k, stacked);
                else
                    imagesc(ax_w, stacked);
                end
                colorbar(ax_w);
                xlabel(ax_w, weight_xlbls{k+1});
                ylabel(ax_w, weight_ylbls{k+1});
                title(ax_w, [part_names{k+1}, ' — weights']);

                % Horizontal separator between components
                hold(ax_w, 'on');
                for c = 1:r_k-1
                    yline(ax_w, c * dim_j + 0.5, 'w-', 'LineWidth', 1);
                end

                % Y-ticks at component midpoints
                mids = (0:r_k-1) * dim_j + round(dim_j / 2) + 1;
                yticks(ax_w, mids);
                yticklabels(ax_w, arrayfun(@(c) sprintf('C%d', c), ...
                    1:r_k, 'UniformOutput', false));

                axis(ax_w, 'tight');  box(ax_w, 'off');
            end
        end

    end % public methods
end % classdef
