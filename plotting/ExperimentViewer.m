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

        %% plotting function headers (normally call external functions that don't rely on custom objects)
        function [hf, out] = plotDistances(obj, ps_lim, method)
            arguments
                obj
                ps_lim double = [1 20] % time interval from stimulus onset [s]
                method char = 'correlation'
            end

            hf = gobjects(1,1);

            % Isolating relevant data
            nsubjects = numel(obj.filtered_traces);
            events = cell(nsubjects,1);
            for i = 1:nsubjects
                thistrace = obj.filtered_traces{i};

                % get peri-stimulus data [t,N,trials]
                M = thistrace.(obj.dataFilter.traceType);
                % # TODO : trial filters
                % # TODO : sorting
                stim_on_frame = thistrace.stim_series.frame_onset(1);
                fs = thistrace.framerate;
                events{i} = TraceViewer.getPeriEventData(M,stim_on_frame,ps_lim,fs);
            end

            % run plotting function
            labs = thistrace.stim_series.stimulus;
            [hf,out] = plotDistances(events,method,labs,obj.plotConfig);
            title(['Similarity: ',method, num2str(ps_lim(1)),'-',num2str(ps_lim(2)), ' s'], ...
                'Color',obj.plotConfig.textcol)

        end

        function obj = setTheme(obj, themeName)
            obj.plotConfig.theme = themeName;
        end
    end
end
