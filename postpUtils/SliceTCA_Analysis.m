classdef SliceTCA_Analysis
% SliceTCA_Analysis  Slice TCA decomposition workflow for neural activity.
%
%   Implements the MATLAB side of a three-step sliceTCA pipeline:
%     1. Export per-subject data tensors to .mat files  (saveInputFiles)
%     2. [External] Run Python sliceTCA script on each file
%     3. [Future]   Load and visualize decomposition     (loadResults)
%
%   The data tensor written per subject follows sliceTCA convention:
%     [trials x neurons x time_bins]  (double)
%
%   Filtering (groups, stimuli, repetitions, time window) is read from
%   obj.v.dataFilter at call time.
%
%   Usage:
%       sla = SliceTCA_Analysis(v);
%       sla.output_dir = 'path/to/slicetca_inputs';
%       paths = sla.saveInputFiles();
%
%       % or pass output_dir directly:
%       paths = sla.saveInputFiles('path/to/slicetca_inputs');

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

    end % public methods
end % classdef
