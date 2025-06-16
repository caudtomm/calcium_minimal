classdef ExperimentViewer
    properties
        name
        subjectTab
        locations
        traces

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
            obj.traces = experiment.traces;
            obj.plotConfig = PlotConfig();
        end

        function out = prepareForPlot(obj)
            % Extract and standardize data for selected subjects and traces

            ids = obj.plotConfig.getSubjectIDs(obj.subjectTab);
            traceTypes = obj.plotConfig.getTraceTypes(obj.traces);

            out = struct( ...
                'subjectIDs', {ids}, ...
                'traceTypes', {traceTypes}, ...
                'traces', [], ...
                'locations', [], ...
                'meta', struct() ...
            );

            numSubjects = numel(ids);
            numTraces = numel(traceTypes);
            out.traces = cell(numTraces, numSubjects);
            out.locations = zeros(numSubjects, size(obj.locations, 2));

            for i = 1:numSubjects
                sid = ids{i};
                idx = find(strcmp({obj.subjectTab.id}, sid), 1);

                if isempty(idx)
                    warning('Subject ID "%s" not found in subjectTab.', sid);
                    continue
                end

                out.locations(i, :) = obj.locations(idx, :);
                for j = 1:numTraces
                    tName = traceTypes{j};
                    if isfield(obj.traces, tName)
                        out.traces{j, i} = obj.traces.(tName)(:, idx);
                    else
                        warning('Trace type "%s" not found in traces.', tName);
                    end
                end
            end

            out.meta.theme = obj.plotConfig.theme;
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
    
                otherwise
                    formatted = raw;
            end
        end

        function obj = setTheme(obj, themeName)
            obj.plotConfig.theme = themeName;
        end
    end
end
