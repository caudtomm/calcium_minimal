function [h,out] = plotTrialMetric(traces,plotType,labs,do_normalize,cfg)
% Low level distance plotter
arguments
    traces cell % cell array [nsubjects 1] of double [1, trials] (sorted!)
    plotType char = 'boxplot'
    labs = [] % stimulus names (sorted!)
    do_normalize logical = false % whether to normalize the traces
    cfg = PlotConfig()
end


% Parse input
assert(numel(labs) == length(traces{1}), 'Mismatch between labs and trial dimension.'); % labs length validation

% Initialize useful metrics
nsubjects = numel(traces);
ntrials = numel(labs);


% Initialize output
out = cell2mat(traces);
h = [];

% Normalize the traces if requested
if do_normalize
    out = out./max(out,[],2,"omitmissing");
end

%% Plot onto provided axes

switch plotType
    case 'boxplot'
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