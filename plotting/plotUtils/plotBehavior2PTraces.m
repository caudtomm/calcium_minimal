function out = plotBehavior2PTraces(tlim, beh_traces, labs, cfg, varargin)
    arguments
        tlim double
        beh_traces cell
        labs cell
        cfg PlotConfig
    end
    arguments (Repeating)
        varargin
    end

    % set default parameters
    dozscorebytrial = false; % zscore each trial individually
    dozscorebysubject = false; % zscore each fish individually

    if ~isempty(varargin)
        for k = 1:2:length(varargin)
            switch lower(varargin{k})
                case 'dozscorebytrial'
                    dozscorebytrial = varargin{k+1};
                case 'dozscorebysubject'
                    dozscorebysubject = varargin{k+1};
            end
        end
    end
    
    % Knobs
    defaultCrange = []; % auto

    if ~isfield(cfg.custom,'crange') || isempty(cfg.custom.crange) % colorbar range as a tunable parameter
        crange = defaultCrange; % fallback
    else
        crange = cfg.custom.crange;
    end

    hasdata = cellfun(@(x) ~isempty(x), beh_traces);
    beh_traces = beh_traces(hasdata);
    labs = labs(hasdata);

    tracesMat = ExperimentViewer.overlayMats(beh_traces);
    [L,nTrials,nSubjects] = size(tracesMat);
    
    y = tracesMat;
    t = linspace(tlim(1), tlim(2), L)';
    labs = labs{1}; % all same labels

    if dozscorebysubject
        for subj = 1:nSubjects
            thismat = y(:,:,subj);
            thismean = mean(thismat, 'all', 'omitmissing');
            thisstd = std(thismat, 0, 'all', 'omitmissing');
            y(:,:,subj) = (thismat - thismean) ./ thisstd;
        end
    end

    if dozscorebytrial; y = nanzscore(y); end
    
    out = mean(y,3,'omitmissing'); % avg over subjects
        
    % plot onto provided axis
    imagesc(1:nTrials,t,out); title('avg over fish'); box off
    xticks(1:nTrials); xticklabels(labs); ylabel('time [s]')
    
    try; clim(crange); catch; end % if crange empty, skip
    colormap(cfg.colormapName)
    a = colorbar('Color',cfg.axcol);
    a.Label.String = 'Value';
    a.Label.FontSize= gca().FontSize;

    set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
    set(gcf, 'color', cfg.bgcol);
end