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
            traceType = obj.plotConfig.traceType;

            out = struct( ...
                'subjectIDs', {ids}, ...
                'traceType', traceType, ...
                'traces', [], ...
                'locations', [], ...
                'meta', struct() ...
            );

            numSubjects = numel(ids);
            out.traces = cell(numSubjects,1);
            out.locations = cell(numSubjects,1);

            for i = 1:numSubjects
                sid = ids{i};
                idx = find(ismember(obj.subjectTab.name, sid), 1);

                if isempty(idx)
                    warning('Subject ID "%s" not found in subjectTab.', sid);
                    continue
                end

                out.locations{i} = obj.traces{idx}.subject_locations;
                out.traces{i} = obj.traces{idx}.(traceType);
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

                case 'plotDistances'
                    % expected input:
                    % .traces - cell array [nsubjects 1] of double [t, cells, trials] (sorted!)
                    % .fs - framerate in seconds (assumes the same for all)
                    % .sec_range - double [absolute absolute] (assumes the same for all)
                    % .labs - cell array [nsubjects 1] of double indices [ntrials 1] with
                    %       stimulus names (sorted!)

                    % # TODO : sorting

                    formatted.traces = raw.traces;
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
