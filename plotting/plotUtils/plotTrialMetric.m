function out = plotTrialMetric(traces,plotType,labs,do_normalize,cfg)
% Low level distance plotter
arguments
    traces cell % cell array [nsubjects 1] of double [1, trials] (sorted!)
    plotType char = 'boxplot'
    labs = [] % stimulus names (sorted!)
    do_normalize logical = false % whether to normalize the traces
    cfg = PlotConfig()
end


% Parse input
assert(numel(labs{1}) == length(traces{1}), 'Mismatch between labs and trial dimension.'); % labs length validation

% Initialize useful metrics
nsubjects = numel(traces);
ntrials = cellfun(@numel,traces);


% Initialize output
out = nan(nsubjects, max(ntrials));
h = [];

% data cell->mat
for i = 1:nsubjects
    n = numel(traces{i});
    out(i, 1:n) = traces{i};
end

% Normalize the traces if requested
if do_normalize
    out = out./max(out,[],2,"omitmissing");
end

%% Plot onto provided axes

switch plotType
    case 'boxplot'
        % select longest labs available
        [~,i]=max(ntrials);
        labs = labs{i}; % #TODO: this is WRONG, because labels may differ

        % plot
        h = boxplot(out, 'Labels', labs, 'PlotStyle', 'compact');
        box off
        set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
        set(gcf, 'color', cfg.bgcol);
    case 'boxplot_repetitions'
        n_repetitions = 5; % # TODO : should be max available repetition`
        out = reshape(out, [], n_repetitions);
        h = boxplot(out, 'PlotStyle', 'compact');
        box off
        set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
        set(gcf, 'color', cfg.bgcol);
    otherwise
        error('Requested plot type is unknown.')
end


end