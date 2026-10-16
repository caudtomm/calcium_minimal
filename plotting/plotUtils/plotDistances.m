function out = plotDistances(traces,plotType,metric,labs,cfg)
% plotDistances - Plots distance matrices for trial data
% 
% Usage:
%   out = plotDistances(traces, plotType, method, labs, cfg)
%
% Inputs:
%   traces   - cell array [nsubjects x 1] of double [t, cells, trials] (sorted)
%   plotType - char, 'full' or 'repetitions'
%   metric   - char/string (pdist method) or DistanceMetric object
%   labs     - cell array, stimulus names (sorted)
%   cfg      - PlotConfig object with plotting parameters
%
% Outputs:
%   out      - distance matrix or similarity matrix depending on plotType
arguments
    traces cell % cell array [nsubjects 1] of double [t, cells, trials] (sorted!)
    plotType char = 'full'
    metric = 'correlation'  % char/string or DistanceMetric
    labs = []               % stimulus names (sorted!)
    cfg = PlotConfig()
end

% Normalize metric argument
if ischar(metric) || isstring(metric)
    metric = DistanceMetric.fromPdist(char(metric));
end

% Knobs
defaultCrange = [.2 1];

% Parse input
assert(numel(labs) == size(traces{1}, 3), 'Mismatch between labs and trial dimension.'); % labs length validation
if ~isfield(cfg.custom,'crange') || isempty(cfg.custom.crange) % colorbar range as a tunable parameter
    crange = defaultCrange; % fallback
else
    crange = cfg.custom.crange;
end

% Initialize useful metrics
nsubjects = numel(traces);
ntrials = numel(labs);


% Initialize output
out=[];

%% Plot onto provided axes

switch plotType
    case 'full' % Plot all in one big matrix
        out = plotFullMat(true);
    case 'repetitions'
        distanceMAT3D = plotFullMat(false);
        out = plotRepetitionSimilarity(...
            distanceMAT3D,unique(labs),true);
    case 'diff_stimulus_repetitions'
        distanceMAT3D = plotFullMat(false);
        out = plotDiffStimSimilarity(...
            distanceMAT3D,unique(labs),true);
    otherwise
        error('Requested plot type is unknown.')
end


%% Functions
function [distanceMAT3D] = plotRepetitionSimilarity(C,stims,pl)
    nstims = numel(stims);
    n_repetitions = 5; % # TODO : should be max available repetition
    distanceMAT3D = nan(n_repetitions,n_repetitions,nsubjects*nstims);
    for i_fish = 1:nsubjects
        c = C(:,:,i_fish);
    
        % now, isolate odors
        for i_stim = 1:nstims
            idx = find(ismember(labs,stims{i_stim})); % only this odor's trials
            distanceMAT3D(:,:,(i_fish-1)*nstims+i_stim) = c(idx,idx);
        end
    end


    % (optional) plot average similarity matrix across repetitions
    if ~pl; return; end
    distanceMAT2D = mean(distanceMAT3D,3,'omitmissing');
    imagesc(1-distanceMAT2D)
    axis square; hold on
    xticks(1:ntrials); xticklabels(1:n_repetitions); xtickangle(90)
    yticks(1:ntrials); yticklabels(1:n_repetitions)
    xlabel('Repetition')
    ylabel('Repetition')
    clim(crange)
    colormap(cfg.colormapName)
    a = colorbar('Color',cfg.axcol);
    a.Label.String = metric.name;
    a.Label.FontSize= gca().FontSize;
    set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
    set(gcf, 'color', cfg.bgcol); 
    hold off
end


function [distanceMAT3D] = plotDiffStimSimilarity(C,stims,pl)
    nstims = numel(stims);
    n_repetitions = 5; % # TODO : should be max available repetition
    if metric.isSymmetric
        n_layers = nstims * (nstims - 1) / 2;
        ks_mask  = ~triu(true(nstims));       % lower triangle only
    else
        n_layers = nstims * (nstims - 1);     % all ordered off-diagonal pairs
        ks_mask  = ~eye(nstims, 'logical');
    end
    distanceMAT3D = nan(n_repetitions, n_repetitions, nsubjects*n_layers);

    k = reshape(1:nstims^2, [nstims, nstims]);
    ks = k(ks_mask);
    k = repelem(k, n_repetitions, n_repetitions);

    n=0;
    for i_fish = 1:nsubjects
        c = C(:,:,i_fish);
    
        for i_k = 1:numel(ks)
            n=n+1;
            distanceMAT3D(:,:,n) = reshape(c(k==ks(i_k)),n_repetitions,n_repetitions);
        end
    end


    % (optional) plot average similarity matrix across repetitions
    if ~pl; return; end
    distanceMAT2D = mean(distanceMAT3D,3,'omitmissing');
    imagesc(1-distanceMAT2D)
    axis square; hold on
    xticks(1:ntrials); xticklabels(1:n_repetitions); xtickangle(90)
    yticks(1:ntrials); yticklabels(1:n_repetitions)
    xlabel('Repetition')
    ylabel('Repetition')
    clim(crange)
    colormap(cfg.colormapName)
    a = colorbar('Color',cfg.axcol);
    a.Label.String = metric.name;
    a.Label.FontSize= gca().FontSize;
    set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
    set(gcf, 'color', cfg.bgcol); 
    hold off
end



function [distanceMAT3D] = plotFullMat(pl)
    distanceMAT3D = nan(ntrials,ntrials,nsubjects);

    for i = 1:nsubjects
        % compress into avg activity vectors
        avg_actvect = squeeze(mean(traces{i},1,'omitmissing')); % [N,trials]
        distanceMAT3D(:,:,i) = metric.compute(avg_actvect'); % [trials, trials]
    end

    % (optional) plot full avg intertrial distance matrix
    if ~pl; return; end
    distanceMAT2D = mean(distanceMAT3D,3,'omitmissing');
    imagesc(1-distanceMAT2D)
    axis square; hold on
    xticks(1:ntrials); xticklabels(labs); xtickangle(90)
    yticks(1:ntrials); yticklabels(labs)
    xlabel('Stimulus type')
    ylabel('Stimulus type')
    clim(crange)
    colormap(cfg.colormapName)
    a = colorbar('Color',cfg.axcol);
    a.Label.String = metric.name;
    a.Label.FontSize= gca().FontSize;
    set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
    set(gcf, 'color', cfg.bgcol); 
    hold off
end

end