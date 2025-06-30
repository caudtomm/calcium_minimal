classdef DataFilter
    properties
       % Selection filters
        subjectIDs cell = {}         % explicit subject ID list (optional)
        subjectGroup = {}            % e.g., a tag or label to select by group;
                                     % options : {'all','familiarized','trained'}
                                     % or {'group1','group2','groupN'}

        % % # TODO : for now, screw all of this. I suspect that individual
        % % function calls are sufficiently idiosyncratic that this actually
        % % makes things more complex, not less.
        % stimuli = {'all_stims'}      % options : {'all_trials','all_stims',
        %                              % 'all CS+','all CS-','all familiar','all novel'} 
        %                              % or {'stim1','stim2','stimN'};
        % time_range double = [0 20]   % time range : [from stim_on, from stim_on] (sec)
        % only_top_variant_units logical = false  % use only top variant units? 
        % order_trials_by char = 'stim_type'      % options: 
        %                                         % {'stim_type' : A-A-A-B-B-B-C-C-C ,
        %                                         % 'trial_num' : A1-A2-B1-A3-B2-B3 ,
        %                                         % 'relative_trial_num' : A1-A2-A3-B1-B2-B3}
        % todo_reltrialnum double = [1:5]

        traceType char = 'dFoverF_good'     % trace type to select (e.g., 'dFoverF_good' or 'pSpike')
                                            % - match name of ActivityTraces property
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

        function ids = getSubjectIDs(obj, subjectTab) 
            % Return list of subject IDs based on filter criteria
            if ~isempty(obj.subjectIDs)
                ids = obj.subjectIDs;
            elseif ~isempty(obj.subjectGroup) && ...
                    ismember('group', subjectTab.Properties.VariableNames)
                allIDs = subjectTab.name;
                match = ismember(subjectTab.group, obj.subjectGroup);
                ids = allIDs(match);
            else
                ids = subjectTab.name;
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