classdef ExperimentViewer
    properties
        name
        subjectTab
        locations
        traces

        subjects_to_use logical
        filtered_traces

        dataFilter DataFilter
        plotConfig PlotConfig
    end

    methods
        function obj = ExperimentViewer(experiment)
            arguments
                experiment Experiment
            end
            obj.name = experiment.name;
            obj.subjectTab = experiment.subjectTab;
            obj.locations = experiment.locations;
            obj.traces = experiment.traces(:);
            obj.dataFilter = DataFilter();
            obj.plotConfig = PlotConfig();
        end

        %% getters

        function idx = get.subjects_to_use(obj)
            subjectIDs = obj.dataFilter.getSubjectIDs(obj.subjectTab);
            idx = ismember(obj.subjectTab.name,subjectIDs);
        end

        function traces = get.filtered_traces(obj)
            % filter subjects
            traces = obj.traces(obj.subjects_to_use);
        end

        %% setters

        function obj = setTheme(obj, themeName)
            obj.plotConfig.theme = themeName;
        end
        
        %% intermediate-level plotting function headers 
        % (normally call external low level functions that don't rely on custom objects)
        % still output a single plot onto provided axes

        function out = plotDistancesHead(obj, varargin)
            % plotDistances - Plot similarity/distance metrics for peri-stimulus traces across subjects.
            %
            % Usage:
            %   [hf, out] = obj.plotDistances('ps_lim', [start end], 'method', methodName, 'trial_sorting', sortingType)
            %
            % Inputs (as name-value pairs):
            %   'ps_lim'        - 2-element vector specifying peri-stimulus window in seconds [default: [1 20]]
            %   'plotType'      - String selecting which plot to produce
            %                     opitons: {'full','repetitions'}
            %   'method'        - String specifying similarity/distance metric (e.g., 'correlation') [default: 'correlation']
            %                     options: see pdist
            %   'trial_sorting' - String specifying trial sorting method (e.g., 'chronological') [default: 'chronological']
            %                     options: {'chronological','random','stim_id'}
            %   'stim_allowed'  - Cell array of strings containing labels of all stimuli to consider. (eg. {'Arg', 'Leu'})
            %                     or stimulus group name [default: 'all trials']
            %                     options: {'all trials','all stimuli','all CS+','all CS-','all familiar','all novel'}
            %
            % Outputs:
            %   out  - Output structure
            arguments
                obj
            end
            arguments (Repeating)
                varargin
            end

            % Set default values
            ps_lim = [1 20];
            plotType = 'full';
            method = 'correlation';
            trial_sorting = 'chronological';
            stim_allowed = 'all trials';

            % Parse name-value pairs
            if ~isempty(varargin)
                for k = 1:2:length(varargin)
                    switch lower(varargin{k})
                        % data filters
                        case 'ps_lim'
                            ps_lim = varargin{k+1};
                        case 'trial_sorting'
                            trial_sorting = varargin{k+1};
                        case 'stims_allowed'
                            stim_allowed = varargin{k+1};
                        % processing parameters
                        case 'method'
                            method = varargin{k+1};
                        % plotting parameters
                        case 'plottype'
                            plotType = varargin{k+1};
                    end
                end
            end

            % initialize output
            out = [];
            
            % Isolating relevant data
            dft = obj.dataFilter;
            dft.interval = ps_lim;
            dft.trial_sorting = trial_sorting;
            dft.repetitions = []; % always use all repetitions
            dft.stims_allowed = stim_allowed;
            [events, all_labs] = dft.filterData(obj);
            if all(cellfun(@isempty,events)); return; end

            % by default, labels are applied based on subject 1
            labs = all_labs{1};

            % Check that all label arrays are the same. This can happen
            % when filtering stimuli according to other properties than
            % identity (ex. trained/novel) from multiple experimental
            % groups.
            all_equal = all(cellfun(@(x) isequal(x, labs), all_labs));
            if ~all_equal
                warning('Mock labels! Actual stimulus identities differ across subjects!')
                % Replace labs with mock labels following the identity structure in labs{1}
                [~, ~, ic] = unique(labs, 'stable');
                mockLabels = arrayfun(@(x) char('A' + x - 1), ic, 'UniformOutput', false);
                labs = mockLabels;
            end

            % call low-level plotter
            out.distMat3d = plotDistances(events,plotType,method,labs,obj.plotConfig);
            title([num2str(ps_lim(1)),'-',num2str(ps_lim(2)), ' s'], ...
                'Color',obj.plotConfig.textcol)

            % return
            out.all_labs = all_labs;

        end


        function out = plotDiscriminationHead(obj, varargin)
            % plotDiscrimination - Plot classification performance / discriminability of a across subjects.
            %
            % Usage:
            %   [hf, out] = obj.plotDistances('ps_lim', [start end], 'method', methodName, 'trial_sorting', sortingType)
            %
            % Inputs (as name-value pairs):
            %   'ps_lim'        - 2-element vector specifying peri-stimulus window in seconds [default: [1 20]]
            %   'plotType'      - String selecting which plot to produce
            %                     opitons: {'full','repetitions'}
            %   'method'        - String specifying similarity/distance metric (e.g., 'correlation') [default: 'correlation']
            %                     options: see pdist
            %   'stim_allowed'  - Cell array of strings containing labels of all stimuli to consider. (eg. {'Arg', 'Leu'})
            %                     or stimulus group name [default: 'all trials']
            %                     options: {'all trials','all stimuli','all CS+','all CS-','all familiar','all novel'}
            %
            % Outputs:
            %   out  - Output structure
            arguments
                obj
            end
            arguments (Repeating)
                varargin
            end

            % Set default values
            ps_lim = [1 20];
            plotType = 'performance_lines';
            method = 'correlation';
            stim_allowed = 'all stimuli';
            reps_touse = [];
            focus_stims = 'all trials'; % by default, no further filtering
            do_zscore = false;

            % Parse name-value pairs
            if ~isempty(varargin)
                for k = 1:2:length(varargin)
                    switch lower(varargin{k})
                        % data filters
                        case 'ps_lim'
                            ps_lim = varargin{k+1};
                        case 'repetitions'
                            reps_touse = varargin{k+1};
                        case 'stims_allowed'
                            stim_allowed = varargin{k+1};
                        % processing parameters
                        case 'method'
                            method = varargin{k+1};
                        % plotting parameters
                        case 'plottype'
                            plotType = varargin{k+1};
                        case 'focus_stims'
                            focus_stims = varargin{k+1};
                        case 'zscore'
                            do_zscore = varargin{k+1};
                    end
                end
            end

            % initialize output
            out = [];
            
            % Isolating relevant data
            dft = obj.dataFilter;
            dft.interval = ps_lim;
            dft.trial_sorting = 'chronological'; % does not matter here
            dft.repetitions = reps_touse;
            dft.stims_allowed = stim_allowed;
            [events, all_labs] = dft.filterData(obj);
            if all(cellfun(@isempty,events)); return; end

            % useful metrics
            nsubjects = numel(obj.filtered_traces);

            % call low-level processor (perform discrimination analysis)
            all_out = cell(nsubjects,1);
            for i = 1:nsubjects
                thisevents = events{i};
                thislabs = all_labs{i};

                % if there are no allowed stimuli here, skip subject
                if isempty(thislabs); continue; end

                % call post-processing function
                all_out{i} = doDiscrimination(thisevents, thislabs, ...
                                                       'method', method);
            end

            % focus on specific trials for plotting (without changing any of the values!)
            focus_trials = cell(nsubjects,1);
            for i = 1:nsubjects
                thistrace = obj.filtered_traces{i};
                thisgroup = thistrace.subject_group;
                desired_stimuli = getStimuliByGroup(thisgroup,focus_stims);
                idx = ismember(all_labs{i}, desired_stimuli);
                focus_trials{i} = find(idx);
            end

            % get rid of empty data
            idx = cellfun(@isempty,all_out) | cellfun(@isempty,focus_trials);
            all_out(idx) = [];
            focus_trials(idx) = [];
            if isempty(all_out); return; end

            % call low-level plotter
            out = plotDiscrimination(all_out,...
                plotType,obj.plotConfig, ...
                'method',method, ...
                'FocusTrials',focus_trials, ...
                'actualreps', reps_touse, ...
                'zscore',do_zscore);

            % return

        end


        function out = plotTrialActivityMetricHead(obj, varargin)
            % plotTrialActivityMetricHead - Plot activity metric (specified as name-value pair argument 'Metric') across subjects.
            %
            % Usage:
            %   [hf, out] = obj.plotActivityMetricHead('ps_lim', [start end], 'method', methodName, 'trial_sorting', sortingType)
            %
            % Inputs (as name-value pairs):
            %   'Metric'        - String specifying activity metric to plot (e.g., 'max-activity', 'dF/F zscore') [default: 'dF/F']
            %   'ps_lim'        - 2-element vector specifying peri-stimulus window in seconds [default: [1 20]]
            %   'plotType'      - String selecting which plot to produce
            %                     opitons: {'full','repetitions'}
            %   'method'        - String specifying similarity/distance metric (e.g., 'correlation') [default: 'correlation']
            %                     options: see pdist
            %   'stim_allowed'  - Cell array of strings containing labels of all stimuli to consider. (eg. {'Arg', 'Leu'})
            %                     or stimulus group name [default: 'all trials']
            %                     options: {'all trials','all stimuli','all CS+','all CS-','all familiar','all novel'}
            %
            % Outputs:
            %   out  - Output structure
            arguments
                obj
            end
            arguments (Repeating)
                varargin
            end

            % Set default values
            ps_lim = [1 20];
            plotType = 'boxplot';
            method = 'participation ratio';
            n_equals = 'cells';
            trial_sorting = 'chronological';
            stim_allowed = 'all stimuli';
            do_normalize = false;

            % Parse name-value pairs
            if ~isempty(varargin)
                for k = 1:2:length(varargin)
                    switch lower(varargin{k})
                        % data filters
                        case 'ps_lim'
                            ps_lim = varargin{k+1};
                        case 'stims_allowed'
                            stim_allowed = varargin{k+1};
                        case 'trial_sorting'
                            trial_sorting = varargin{k+1};
                        % processing parameters
                        case 'method'
                            method = varargin{k+1};
                        case 'n_equals'
                            n_equals = varargin{k+1};
                        % plotting parameters
                        case 'plottype'
                            plotType = varargin{k+1};
                        case 'normalize'
                            do_normalize = varargin{k+1};
                    end
                end
            end

            % initialize output
            out = [];
            
            % Isolating relevant data
            dft = obj.dataFilter;
            dft.interval = ps_lim;
            dft.trial_sorting = trial_sorting;
            dft.repetitions = []; % always use all repetitions
            dft.stims_allowed = stim_allowed;
            [events, all_labs] = dft.filterData(obj);
            if all(cellfun(@isempty,events)); return; end

            % useful metrics
            nsubjects = numel(obj.filtered_traces);

            % call low-level processor (perform discrimination analysis)
            all_out = cell(nsubjects,1);
            for i = 1:nsubjects
                thisevents = events{i};

                % call post-processing function
                all_out{i} = extractActivityMetric(thisevents, method, n_equals);
            end

            % get rid of empty data
            idx = cellfun(@isempty,all_out);
            all_out(idx) = [];
            if isempty(all_out); return; end

            % if output metrics are not 1D, assume that the last available dimension is trials
            for i = 1:numel(all_out)
                sz = size(all_out{i});
                if numel(sz) > 1 && sz(end) > 1
                    % average along all other dimensions except the last (trials)
                    dims_to_avg = 1:(numel(sz)-1);
                    all_out{i} = squeeze(mean(all_out{i}, dims_to_avg));
                    % ensure column vector
                    if isrow(all_out{i})
                        all_out{i} = all_out{i}';
                    end
                end
            end

            % call low-level plotter
            out = plotTrialMetric(all_out,...
                plotType, all_labs, ...
                do_normalize, obj.plotConfig);

            % return

        end


        function out = plotUnitActivityMetricHead(obj, varargin)
            % plotUnitActivityMetricHead - Plot classification performance / discriminability of a across subjects.
            %
            % Usage:
            %   [hf, out] = obj.plotDistances('ps_lim', [start end], 'method', methodName, 'trial_sorting', sortingType)
            %
            % Inputs (as name-value pairs):
            %   'ps_lim'        - 2-element vector specifying peri-stimulus window in seconds [default: [1 20]]
            %   'plotType'      - String selecting which plot to produce
            %                     opitons: {'full','repetitions'}
            %   'method'        - String specifying similarity/distance metric (e.g., 'correlation') [default: 'correlation']
            %                     options: see pdist
            %   'stim_allowed'  - Cell array of strings containing labels of all stimuli to consider. (eg. {'Arg', 'Leu'})
            %                     or stimulus group name [default: 'all trials']
            %                     options: {'all trials','all stimuli','all CS+','all CS-','all familiar','all novel'}
            %
            % Outputs:
            %   out  - Output structure
            arguments
                obj
            end
            arguments (Repeating)
                varargin
            end

            % Set default values
            ps_lim = [1 20];
            plotType = 'boxplot';
            method = 'tuning curves';
            n_equals = 'cells';
            trial_sorting = 'chronological';
            stim_allowed = 'all stimuli';
            do_normalize = false;

            % Parse name-value pairs
            if ~isempty(varargin)
                for k = 1:2:length(varargin)
                    switch lower(varargin{k})
                        % data filters
                        case 'ps_lim'
                            ps_lim = varargin{k+1};
                        case 'stims_allowed'
                            stim_allowed = varargin{k+1};
                        case 'trial_sorting'
                            trial_sorting = varargin{k+1};
                        % processing parameters
                        case 'method'
                            method = varargin{k+1};
                        case 'n_equals'
                            n_equals = varargin{k+1};
                        % plotting parameters
                        case 'plottype'
                            plotType = varargin{k+1};
                        case 'normalize'
                            do_normalize = varargin{k+1};
                    end
                end
            end

            % initialize output
            out = [];
            
            % Isolating relevant data
            dft = obj.dataFilter;
            dft.interval = ps_lim;
            dft.trial_sorting = trial_sorting;
            dft.repetitions = []; % always use all repetitions
            dft.stims_allowed = stim_allowed;
            [events, all_labs] = dft.filterData(obj);
            if all(cellfun(@isempty,events)); return; end

            % useful metrics
            nsubjects = numel(obj.filtered_traces);

            % call low-level processor (perform discrimination analysis)
            all_out = cell(nsubjects,1);
            for i = 1:nsubjects
                thisevents = events{i};
                thislabs = all_labs{i};

                % if there are no allowed stimuli here, skip subject
                if isempty(thislabs); continue; end

                % call post-processing function
                all_out{i} = extractActivityMetric(thisevents, method, n_equals, ...
                    'StimTypes', thislabs);
            end

            % get rid of empty data
            idx = cellfun(@isempty,all_out);
            all_out(idx) = [];
            if isempty(all_out); return; end

            % call low-level plotter
            % out = plotUnitActivityMetric(all_out,...
            %     plotType, all_labs, ...
            %     do_normalize, obj.plotConfig);

            % return
            out = all_out;

        end



        %% complex and idiosyncratic high-level plotters are saved in external files


    end
end
