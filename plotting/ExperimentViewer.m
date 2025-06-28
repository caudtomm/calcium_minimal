classdef ExperimentViewer
    properties
        name
        subjectTab
        locations
        traces

        subjects_to_use logical

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

        %% functional methods

        function out = prepareForPlot(obj)
            % Extract and standardize data for selected subjects and traces

            ids = obj.subjectTab.name(obj.subjects_to_use);
            traceType = obj.dataFilter.traceType;

            idx = obj.subjects_to_use;

            out = struct( ...
                'subjectIDs', ids, ...
                'traceType', traceType, ...
                'traces', cellfun(@(x) x.(traceType),obj.traces(idx),'UniformOutput',false) ...
            );
        end

        function plotSummary(obj, plotFunc)
            % Call an external plot function using formatted summary data
            raw = obj.prepareForPlot();
            formatted = obj.formatFor(plotFunc, raw);
            plotFunc(formatted);
        end

        function formatted = formatFor(obj, plotFunc, raw)
            fname = func2str(plotFunc);

            switch fname
                case 'plotMeanTrace' % just an example
                    % Assume single traceType is selected
                    traces = raw.traces(1, :);  % row: traceType=1, all subjects
                    valid = ~cellfun(@isempty, traces);
                    dataMat = cell2mat(traces(valid)');
                    meanTrace = mean(dataMat, 2);
    
                    formatted.meanTrace = meanTrace;
                    formatted.subjectIDs = raw.subjectIDs(valid);
                    formatted.traceType = raw.traceTypes{1};
                    formatted.themeColors = struct( ...
                        'axis', obj.plotConfig.getAxisColor(), ...
                        'text', obj.plotConfig.getTextColor(), ...
                        'background', obj.plotConfig.getBackgroundColor(), ...
                        'line', obj.plotConfig.getColorCycle(1) ...
                    );

                case 'plotDistances'
                    % expected input:
                    % .traces - cell array [nsubjects 1] of double [t, cells, trials] (sorted!)
                    % .fs - framerate in seconds (assumes the same for all)
                    % .sec_range - double [absolute absolute] (assumes the same for all)
                    % .labs - cell array [nsubjects 1] of double indices [ntrials 1] with
                    %       stimulus names (sorted!)

                    % # TODO : sorting

                    formatted.traces = raw.traces;

                    % # TODO : individual assignment
                    formatted.fs = obj.traces{1}.framerate;
                    formatted.sec_range = [30 50];
                    formatted.labs = obj.traces{1}.stim_series.stimulus;
                    
    
                otherwise
                    formatted = raw;
            end
        end

        function obj = setTheme(obj, themeName)
            obj.plotConfig.theme = themeName;
        end
    end
end
