function out = plotDistances(traces,plotType,method,labs,cfg)
% Low level distance plotter
arguments
    traces cell % cell array [nsubjects 1] of double [t, cells, trials] (sorted!)
    plotType char = 'full'
    method char = 'correlation'
    labs = [] % stimulus names (sorted!)
    cfg = PlotConfig()
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
    a.Label.String = method;
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
        distanceMAT3D(:,:,i) = squareform(pdist(avg_actvect', method)); % [trials, trials]
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
    a.Label.String = method;
    a.Label.FontSize= gca().FontSize;
    set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
    set(gcf, 'color', cfg.bgcol); 
    hold off
end

end