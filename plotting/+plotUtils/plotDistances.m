function [hf, out] = plotDistances(traces,method,labs,cfg)
arguments
    traces cell % cell array [nsubjects 1] of double [t, cells, trials] (sorted!)
    method char = 'correlation'
    labs = [] % stimulus names (sorted!)
    cfg = PlotConfig()
end

% Knobs
crange = [.2 1];


% Initialize useful metrics
nsubjects = numel(traces);
ntrials = numel(labs);


% Initialize output
hf = gobjects(1,1);
out.distanceMAT3D=[];


% Plot all in one big matrix
[hf(1), out.distanceMAT3D] = plotFullMat(true);






function [hf, distanceMAT3D] = plotFullMat(pl)
    distanceMAT3D = nan(ntrials,ntrials,nsubjects);

    for i = 1:nsubjects
        % compress into avg activity vectors
        avg_actvect = squeeze(mean(traces{i},1,'omitmissing')); % [N,trials]
        distanceMAT3D(:,:,i) = squareform(pdist(avg_actvect', method)); % [trials, trials]
    end

    % (optional) plot full avg intertrial distance matrix
    if ~pl; return; end
    distanceMAT2D = mean(distanceMAT3D,3,'omitmissing');
    hf = figure('Color',cfg.bgcol); imagesc(1-distanceMAT2D)
    axis square; hold on
    xticks(1:ntrials); xticklabels(labs); xtickangle(90)
    yticks(1:ntrials); yticklabels(labs)
    xlabel('Stimulus type')
    ylabel('Stimulus type')
    clim(crange)
    colormap(cfg.colormapName)
    colorbar('Color',cfg.axcol)
    set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
    set(gcf, 'color', cfg.bgcol); 
    hold off
end

end