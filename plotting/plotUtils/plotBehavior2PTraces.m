function out = plotBehavior2PTraces(t, beh_traces, labs, cfg, varargin)
    arguments
        t cell
        beh_traces cell
        labs cell
        cfg PlotConfig
    end
    arguments (Repeating)
        varargin
    end

    hasdata = cellfun(@(x) ~isempty(x), beh_traces);
    beh_traces = beh_traces(hasdata);
    t = t(hasdata);
    labs = labs(hasdata);

    tracesMat = ExperimentViewer.overlayMats(beh_traces);
    
    y = tracesMat; % time x labels x subjects
    t = t{1}; % all same time axis
    labs = labs{1}; % all same labels

    out = mean(nanzscore(y),3,'omitmissing'); % avg over subjects
        
    % plot onto provided axis
    imagesc(1:35,t,out); title('avg over fish'); box off
    xticks(1:35); xticklabels(labs); ylabel('time [s]')
    set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
    set(gcf, 'color', cfg.bgcol);
end