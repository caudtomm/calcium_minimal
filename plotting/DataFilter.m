classdef DataFilter
    properties
       % Selection filters
        subjectIDs cell = {}         % explicit subject ID list (optional)
        subjectGroup = {}            % e.g., a tag or label to select by group;
                                     % options : {'all','familiarized','trained'}
                                     % or {'group1','group2','groupN'}

        traceType char = 'dFoverF_good'     % trace type to select (e.g., 'dFoverF_good' or 'pSpike')
                                            % - match name of ActivityTraces property

        interval double = [1 20] % peristimulus limits [sec]
        trial_sorting char = 'stim_id' % trial sorting method, options: {'stim_id', 'chronological', 'relative_trial_num'}
        repetitions double = [] % stimulus repetitions to use (empty = all)
        stims_allowed = 'all stimuli' % list of allowed stimuli, type cell or char vector, see accepted inputs to getStimuliByGroup()

        % behavior2p related properties
        behavior2p_trace char = 'tail_motion_2p' % behavior2p trace type to select{'breathing_events', 'breathing_ipis', ...
                            % 'breathing_inst_freq', 'tail_motion', 'breathing_inst_freq_2p', 'tail_motion_2p'}
        
        % mode selection (unused here, but passed on to ModeSelector)
        mode_name char = 'native_units' % {'native_units', 'pca', 'nmf', 'dpca'}
        mode_method char = 'mode_values' % {'mode_values', 'isolate', 'subtract'}
        mode_OI = 'all' % modes of interest: {'all', 'stimulus', 'non-stimulus', 'novelty'}
        mode_file char = ''; % if empty, extracts by default. Else, it looks for coefficients in the file specified.
        mode_params = struct(); % specify manually as needed
    end

    methods (Static)
        function groups_to_use = parseGroupTag(group)
            % Static method to parse a group tag string
            % Group names can be char vectors, strings or cell arrays of
            % strings;
            
            % Default group sets tags:
            defaultSets = {'all', 'familiarized', 'trained'};

            % List of group names found in the dataset
            knownGroups = {'previousnaive'; ...
                           'naïve'; ...
                           'trained1'; ...
                           'trained2'; ...
                           'trained1 (T-R-S-H-A-ACSF/L)'; ...
                           'uncoupled'};

            % Parse group tag string to a standardized format
            if ischar(group) || isstring(group) % single group name or group set
                group = {char(group)};
            elseif ~iscell(group)
                error('Group tag must be a string or cell array of strings.');
            end

            % Now group is cell
            
            % Check valid names
            validGroups = cellfun(@isValidGroup, group);
            if ~all(validGroups)
                missing = group(~validGroups);
                warning('Unknown group(s): %s', strjoin(missing, ', '));
            end
            group = group(validGroups);

            % Go through list
            tmp = cellfun(@expandSet, group, 'UniformOutput', false);
            groups_list = vertcat(tmp{:});

            % Return result
            groups_to_use = unique(groups_list);

            % Functions
            function isValid = isValidGroup(g)
                % Check if group is a valid single group or a default set
                isValid = any(strcmp(g, knownGroups)) || ...
                          any(strcmp(g, defaultSets));
            end

            function g_list = expandSet(s)
                switch s
                    case 'all'
                        g_list = knownGroups; % all known groups
                    case 'familiarized'
                        % subjects who have been previously exposed to some of the stimuli
                        g_list = {'trained1'; ...
                                 'trained2'; ...
                                 'trained1 (T-R-S-H-A-ACSF/L)'; ...
                                 'uncoupled'};
                    case 'trained'
                        % subjects who have been trained on the task
                        g_list = {'trained1'; ...
                                 'trained2'; ...
                                 'trained1 (T-R-S-H-A-ACSF/L)'};
                    otherwise   % must be one of the single group names:
                        % can be assumed because we checked for validity
                        % earlier.
                        g_list = {s};
                end
            end
        end
    end

    methods
        function obj = DataFilter(varargin)        
            if nargin == 1 && isstruct(varargin{1})
                % Initialize from struct
                s = varargin{1};
                fn = fieldnames(s);
                for i = 1:numel(fn)
                    if isprop(obj, fn{i})
                        obj.(fn{i}) = s.(fn{i});
                    end
                end
            elseif mod(nargin, 2) == 0
                % Name-value pair input
                for i = 1:2:nargin
                    name = varargin{i};
                    value = varargin{i+1};
                    if isprop(obj, name)
                        obj.(name) = value;
                    else
                        error('Invalid property name: %s', name);
                    end
                end
            elseif nargin > 0
                error('Unsupported DataFilter constructor usage.');
            end
        end

        %% getters

        function [ids, match] = getSubjectIDs(obj, subjectTab) 
            % Return list of subject IDs based on filter criteria
            if ~isempty(obj.subjectIDs)
                ids = obj.subjectIDs;
                match = true(numel(ids),1);
            elseif ~isempty(obj.subjectGroup) && ...
                    ismember('group', subjectTab.Properties.VariableNames)
                allIDs = subjectTab.name;
                match = ismember(subjectTab.group, obj.subjectGroup);
                ids = allIDs(match);
            else
                ids = subjectTab.name;
                match = true(numel(ids),1);
            end
        end


        %% setters
        function obj = set.subjectIDs(obj, ids)
            if iscell(ids) && all(cellfun(@ischar, ids))
                obj.subjectIDs = ids(:);
            else
                error('subjectIDs must be a cell array of strings.');
            end
        end

        function obj = set.subjectGroup(obj, group)
            obj.subjectGroup = obj.parseGroupTag(group);
        end

        %% filtering

        function [events, labs] = filterData(obj, v)
            arguments
                obj (1,1) DataFilter
                v (1,1) ExperimentViewer
            end
            
            % extract properties to vars
            traceType = obj.traceType; % trace type to select (e.g., 'dFoverF_good' or 'pSpike')
            ps_lim = obj.interval; % peristimulus limits [sec]
            reps_touse = obj.repetitions;
            stim_allowed = obj.stims_allowed;
            trial_sorting = obj.trial_sorting;

            % Filter data based on the properties of this DataFilter object
            nsubjects = numel(v.filtered_traces);
            events = cell(nsubjects,1);
            labs = cell(nsubjects,1);
            for i = 1:nsubjects
                thistrace = v.filtered_traces{i};

                % Trial sorting
                [~,trial_idx] = TraceViewer(thistrace).sortTrials(trial_sorting);

                % get peri-stimulus data [t,N,trials]
                M = thistrace.(traceType)(:,:,trial_idx);
                stim_on_frame = thistrace.stim_series.frame_onset(1);
                fs = thistrace.framerate;
                events{i} = TraceViewer.getPeriEventData(M,stim_on_frame,ps_lim,fs);

                % Retrieve stimulus identity labels
                labs{i} = thistrace.stim_series.stimulus(trial_idx);

                % Stimulus filtering
                thisgroup = thistrace.subject_group;
                desired_stimuli = getStimuliByGroup(thisgroup,stim_allowed);
                idx = ismember(labs{i}, desired_stimuli);
                if ~isempty(desired_stimuli) && ~all(idx)
                    % If some trials are not in the desired stimuli, filter them out
                    events{i} = events{i}(:,:,idx);
                    labs{i} = labs{i}(idx);
                elseif isempty(desired_stimuli) || sum(idx)==0
                    % If no stimuli are accepted, return empty arrays!
                    events{i} = [];
                    labs{i} = [];
                    return
                end

                % Stimulus repetition filter
                if isempty(reps_touse); continue; end % empty argument 'repetitions' leads to all repetitions being used
                thisstims = unique(labs{i});
                nstims = numel(thisstims);
                idx_keep = false(1, numel(labs{i}));
                for i_stim = 1:nstims
                    idx_stim = find(ismember(labs{i}, thisstims{i_stim}));
                    this_nreps = numel(idx_stim);
                    % Select only allowed repetition indices
                    reps_available = 1:this_nreps;
                    reps_valid = reps_available(ismember(reps_available, reps_touse));
                    if isempty(reps_valid)
                        continue
                    end
                    idx_keep(idx_stim(reps_valid)) = true;
                end
                if ~all(idx_keep)
                    % Filter events and labels to keep only desired repetitions
                    events{i} = events{i}(:,:,idx_keep);
                    labs{i} = labs{i}(idx_keep);
                elseif sum(idx_keep)==0
                    % If no trials are accepted, return without trying to
                    % plot ... nothing!
                    return
                end
            end

        end

        %% export
        function s = export(obj)
            % Optional: dump config as a struct
            s = struct( ...
                'subjectIDs', obj.subjectIDs, ...
                'subjectGroup', obj.subjectGroup, ...
                'traceType', obj.traceType ...
            );
        end
    end

end