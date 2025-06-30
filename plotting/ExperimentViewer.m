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

            % % filtering within each subject traces
            % nsubjects = numel(traces);
            % for i = 1:nsubjects
            %     % filter 
            % end
        end

        %% setters

        function obj = setTheme(obj, themeName)
            obj.plotConfig.theme = themeName;
        end
        
        %% intermediate-level plotting function headers 
        % (normally call external low level functions that don't rely on custom objects)
        % still output a single plot onto provided axes

        function out = plotDistances(obj, varargin)
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
                        case 'ps_lim'
                            ps_lim = varargin{k+1};
                        case 'plottype'
                            plotType = varargin{k+1};
                        case 'method'
                            method = varargin{k+1};
                        case 'trial_sorting'
                            trial_sorting = varargin{k+1};
                        case 'stim_allowed'
                            stim_allowed = varargin{k+1};
                    end
                end
            end

            % initialize output
            out = [];

            % Isolating relevant data
            nsubjects = numel(obj.filtered_traces);
            events = cell(nsubjects,1);
            all_labs = cell(nsubjects,1);
            for i = 1:nsubjects
                thistrace = obj.filtered_traces{i};

                % Trial sorting
                [~,trial_idx] = TraceViewer(thistrace).sortTrials(trial_sorting);

                % get peri-stimulus data [t,N,trials]
                M = thistrace.(obj.dataFilter.traceType)(:,:,trial_idx);
                stim_on_frame = thistrace.stim_series.frame_onset(1);
                fs = thistrace.framerate;
                events{i} = TraceViewer.getPeriEventData(M,stim_on_frame,ps_lim,fs);

                % Retrieve stimulus identity labels
                all_labs{i} = thistrace.stim_series.stimulus(trial_idx);

                % Stimulus filtering
                thisgroup = thistrace.subject_group;
                desired_stimuli = getStimuliByGroup(thisgroup,stim_allowed);
                idx = ismember(all_labs{i}, desired_stimuli);
                if ~isempty(desired_stimuli) && ~all(idx) && sum(idx)>0
                    % If some trials are not in the desired stimuli, filter them out
                    events{i} = events{i}(:,:,idx);
                    all_labs{i} = all_labs{i}(idx);
                elseif isempty(desired_stimuli) || sum(idx)==0
                    % If no stimuli are accepted, return without trying to
                    % plot ... nothing!
                    return
                end
            end

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
            out.distMat3d = plotUtils.plotDistances(events,plotType,method,labs,obj.plotConfig);
            title([num2str(ps_lim(1)),'-',num2str(ps_lim(2)), ' s'], ...
                'Color',obj.plotConfig.textcol)

            % return
            out.all_labs = all_labs;

        end


        function out = plotDiscrimination(obj, varargin)
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
                        case 'ps_lim'
                            ps_lim = varargin{k+1};
                        case 'plottype'
                            plotType = varargin{k+1};
                        case 'method'
                            method = varargin{k+1};
                        case 'trial_sorting'
                            trial_sorting = varargin{k+1};
                        case 'stim_allowed'
                            stim_allowed = varargin{k+1};
                    end
                end
            end

            % initialize output
            out = [];

            % Isolating relevant data
            nsubjects = numel(obj.filtered_traces);
            events = cell(nsubjects,1);
            all_labs = cell(nsubjects,1);
            for i = 1:nsubjects
                thistrace = obj.filtered_traces{i};

                % Trial sorting
                [~,trial_idx] = TraceViewer(thistrace).sortTrials(trial_sorting);

                % get peri-stimulus data [t,N,trials]
                M = thistrace.(obj.dataFilter.traceType)(:,:,trial_idx);
                stim_on_frame = thistrace.stim_series.frame_onset(1);
                fs = thistrace.framerate;
                events{i} = TraceViewer.getPeriEventData(M,stim_on_frame,ps_lim,fs);

                % Retrieve stimulus identity labels
                all_labs{i} = thistrace.stim_series.stimulus(trial_idx);

                % Stimulus filtering
                thisgroup = thistrace.subject_group;
                desired_stimuli = getStimuliByGroup(thisgroup,stim_allowed);
                idx = ismember(all_labs{i}, desired_stimuli);
                if ~isempty(desired_stimuli) && ~all(idx)
                    % If some trials are not in the desired stimuli, filter them out
                    events{i} = events{i}(:,:,idx);
                    all_labs{i} = all_labs{i}(idx);
                elseif isempty(desired_stimuli) || sum(idx)==0
                    % If no stimuli are accepted, return without trying to
                    % plot ... nothing!
                    return
                end
            end

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
            out.distMat3d = plotUtils.plotDistances(events,plotType,method,labs,obj.plotConfig);
            title([num2str(ps_lim(1)),'-',num2str(ps_lim(2)), ' s'], ...
                'Color',obj.plotConfig.textcol)

            % return
            out.all_labs = all_labs;

        end



        %% complex and idiosyncratic high-level plotters

        function [hf, out] = plotRepetitionDistances(obj, ps_lim, method)  
            arguments
                obj
                ps_lim = [1 20]
                method = 'correlation'
            end
            
            % Knobs
            groups = {'naïve', ...
                      'trained', ...
                      'trained1', ...
                      'trained2', ...
                      'trained1 (T-R-S-H-A-ACSF/L)', ...
                      'uncoupled'};
            odor_sets = {'all stimuli', ...
                         {'Arg','Ala','His','Trp','Ser'}, ...
                         {'Leu'}, ...
                         'all CS+', ...
                         'all CS-', ...
                         'all familiar', ...
                         'all novel', ...
                        };

            % useful metrics
            ngroups = numel(groups);
            nodor_sets = numel(odor_sets);
            nplots = ngroups*nodor_sets;

            % Initialize output
            hf = gobjects(2,1);
            out = cell(nplots,1);


            %% Figure 1: comparison of distance matrices

            % define figure size
            ncols = nodor_sets;
            nrows = ngroups;

            hf(1) = figure; % [groups, odor_sets]
            n = 1;
            labs = cell(nplots,1);
            for i_g = 1:ngroups
                % filter data by group
                thisgroup = groups{i_g};
                obj.dataFilter.subjectGroup = thisgroup;

                for i_o = 1:nodor_sets
                    thisodorset = odor_sets{i_o};

                    thisodorset_str = thisodorset;
                    if iscell(thisodorset_str); thisodorset_str = strjoin(thisodorset, ', '); end
                    msg = ['Plotting group ''',thisgroup,''' for odors: ',thisodorset_str];
                    disp(msg)

                    % build axes
                    subplot(nrows,ncols,n);

                    % call intermediate-level plotter
                    out{n} = obj.plotDistances('ps_lim',ps_lim, ...
                                    'plotType', 'repetitions', ...
                                    'method',method, ...
                                    'stim_allowed',thisodorset);

                    % override title and ylabel
                    title(thisodorset_str)
                    ylabel(thisgroup)

                    % if there is no data, delete the subplot
                    if isempty(out{n}); axis off; end

                    % export label
                    labs{n} = [thisgroup,' - ',thisodorset_str];

                    % advance axis counter
                    n = n+1;
                end
            end
            set(gcf,'Position',[1 1 2000 1000])

            %% Figure 2: comparison of distance distributions

            % compute data
            data = [];
            for i = 1:nplots
                % skip empty output
                if isempty(out{i}); continue; end

                thismat = 1 - out{i}.distMat3d; % convert distance to similarity
                
                % knobs # TODO : tunable param
                repetitions = 2:5;

                % crop and take out the diagonal
                thismat = thismat(repetitions,repetitions,:);
                identity = logical(eye(length(repetitions)));
                thismat(identity) = nan;

                % store
                plot_idx = i * ones(numel(thismat),1);
                data = [data; plot_idx, thismat(:)];
            end
            
            % plot
            cfg = obj.plotConfig;
            hf(2) = figure;
            boxplot(data(:,2),labs(data(:,1)),'Orientation','horizontal','PlotStyle','compact','Colors',cfg.textcol)
            box off
            set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol);
            set(gcf, 'color', cfg.bgcol); 
            xlabel(method,'Color',cfg.textcol)
            set(gcf,'Position',[1 1 2000 1000])
            hold off
        end
    end
end
