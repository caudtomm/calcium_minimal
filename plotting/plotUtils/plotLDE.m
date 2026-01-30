function out = plotLDE(embedding,plotType,labs,cfg, varargin)
% <Low level low dimensional embedding plotter>
arguments
    embedding double % [nSamples x nDimensions x nTrials] matrix
    plotType char = 'scatter'
    labs = [] % stimulus names (sorted!)
    cfg = PlotConfig()
end
arguments (Repeating)
    varargin
end

% Default parameters
params.ldeType = 'projection';
params.sizeByRep = true; % size points by repetition number

% Parse name-value pairs
if ~isempty(varargin)
    for k = 1:2:length(varargin)
        switch lower(varargin{k})
            case 'ldetype'
                params.ldeType = varargin{k+1};
            otherwise
                error('Unknown parameter name: %s', varargin{k});
        end
    end
end

% Parse input
assert(numel(labs) == size(embedding, 3), 'Mismatch between labs and trial dimension.'); % labs length validation

% Initialize output
out=[];

%% Plot onto provided axes
switch plotType
    case 'scatter' % Plot all in one big matrix
        out = plotScatter(embedding,labs,params,cfg);
    case 'lines'
        out = plotLines(embedding,labs,params,cfg);
    otherwise
        error('Requested plot type is unknown.')
end

end


%% Functions
function h = plotScatter(embedding,labs,params,cfg)
    [L, ndims, ntrials] = size(embedding);

    % set scatter size
    scalefact = 500;
    sizesc = 10+scalefact/(L*ntrials);

    stims = unique(labs);

    h = gobjects(1, ntrials); % Preallocate h array for better performance
    for i_trial = 1:ntrials
        thisStim = labs(i_trial);
        thisStimID = find(ismember(stims,thisStim));
        thisRepNum = sum(ismember(labs(1:i_trial),thisStim));

        if numel(stims) > 1
            thisColor = cfg.c(thisStimID,:); % color by stimulus
        else
            thisColor = cfg.c(thisRepNum,:); %color by repetition
        end

        sizeMultiplier = 1;
        if params.sizeByRep && numel(stims) > 1
            sizeMultiplier = 1/thisRepNum;
        end

        switch ndims
            case 2
                h(i_trial) = scatter(embedding(:,1,i_trial),embedding(:,2,i_trial), ...
                    sizesc*sizeMultiplier,thisColor,'filled');
            case 3
                h(i_trial) = scatter3(embedding(:,1,i_trial),embedding(:,2,i_trial), ...
                    embedding(:,3,i_trial),sizesc*sizeMultiplier, ...
                    thisColor,'filled');
            otherwise
                error('dimension number not supported')
        end
        
        if numel(stims) > 1
            h(i_trial).DisplayName = thisStim{1}; % legend by stimulus
        else
            h(i_trial).DisplayName = ['Rep #',num2str(thisRepNum)]; % legend by repetition
        end
        hold on
    end

    
    % cosmetics / labels
    u = legendUnq();
    legend(u,'Box','on','color',cfg.bgcol,'Location','best', ...
        'EdgeColor',cfg.textcol,'TextColor',cfg.textcol)
    set(gcf, 'color', cfg.bgcol);    
    set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, ...
        'YColor',cfg.axcol, 'ZColor',cfg.axcol);
    xlabel([params.ldeType,' #1'])
    ylabel([params.ldeType,' #2'])
    if ndims>2; zlabel([params.ldeType,' #3']); end
    grid off
    axis square
end

function h = plotLines(embedding,labs,params,cfg)
    [L, ndims, ntrials] = size(embedding);

    % set scatter size
    scalefact = 5;
    sizesc = 10+scalefact/(L*ntrials);

    stims = unique(labs);

    for i_trial = 1:ntrials
        thisStim = labs(i_trial);
        thisStimID = find(ismember(stims,thisStim));
        thisRepNum = sum(ismember(labs(1:i_trial),thisStim));

        thisColor = cfg.c(thisStimID,:);

        sizeMultiplier = 1;
        if params.sizeByRep
            sizeMultiplier = .3/thisRepNum;
        end

        switch ndims
            case 2
                h(i_trial) = plot(embedding(:,1,i_trial),embedding(:,2,i_trial), ...
                    'LineWidth',sizesc*sizeMultiplier,'Color',thisColor);
            case 3
                h(i_trial) = plot3(embedding(:,1,i_trial),embedding(:,2,i_trial), ...
                    embedding(:,3,i_trial), ...
                    'LineWidth',sizesc*sizeMultiplier, ...
                    'Color',thisColor);
            otherwise
                error('dimension number not supported')
        end
        h(i_trial).DisplayName = thisStim{1};
        hold on
    end

    
    % cosmetics / labels
    u = legendUnq();
    legend(u,'Box','on','color',cfg.bgcol,'Location','best', ...
        'EdgeColor',cfg.textcol,'TextColor',cfg.textcol)
    set(gcf, 'color', cfg.bgcol);    
    set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, ...
        'YColor',cfg.axcol, 'ZColor',cfg.axcol);
    xlabel([params.ldeType,' #1'])
    ylabel([params.ldeType,' #2'])
    if ndims>2; zlabel([params.ldeType,' #3']); end
    grid off
    axis square
end
