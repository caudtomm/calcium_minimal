function hf = compareModeMetricsFigure(v, dft)
    cfg = v.plotConfig;
    
    m = ModeSelector(v,'native_units',dft); % extract modes
    [GSS, SeT, StT, Tuning] = getSelectedMetrics(m);
    hf = ScatterFigure();
    
    % pick 5 random subjects
    nsubjects = numel(m.data);
    nsubjects2use = 5;
    idx = randi(nsubjects,nsubjects2use,1);
    
    hf = ExamplesFigure(idx);
    
    % NMF
    m = ModeSelector(v,'nmf',dft,'nfactors',15); % extract modes
    [GSS, SeT, StT, Tuning] = getSelectedMetrics(m);
    hf = ScatterFigure();
    hf = ExamplesFigure(idx);
    
    % PCA
    m = ModeSelector(v,'pca',dft,'nfactors',15); % extract modes
    [GSS, SeT, StT, Tuning] = getSelectedMetrics(m);
    hf = ScatterFigure();
    hf = ExamplesFigure(idx);

    %% functions

    function [GSS, SeT, StT, Tuning] = getSelectedMetrics(m)
        labels = m.labels;
        events = m.values;
        GSS = loopMetricExtraction(events,labels,'general suppression score');
        SeT = loopMetricExtraction(events,labels,'selectivity of tuning');
        StT = loopMetricExtraction(events,labels,'stability of tuning');
        
        Tuning = loopMetricExtraction(events,labels,'tuning curves');
        
        % GSS = cell2mat(GSS'); % [modes x subjects]
        % SeT = cell2mat(SeT'); % [modes x subjects]
        % StT = cell2mat(StT'); % [modes x subjects]
    end
    
    function hf = ScatterFigure()
        hf = figure;
        
        subplot(131)
        plotScatter(StT,GSS,cfg,'subject')
        box off
        axis square
        xlabel('Stability of Tuning')
        ylabel('g. Suppression')
        
        subplot(132)
        plotScatter(SeT,GSS,cfg,'subject')
        box off
        axis square
        xlabel('Selectivity of Tuning')
        ylabel('g. Suppression')
        
        subplot(133)
        plotScatter(StT,SeT,cfg,'subject')
        box off
        axis square
        xlabel('Stability of Tuning')
        ylabel('Selectivity of Tuning')
        
        for i = 1:3
            subplot(1,3,i)
            set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol);
        end
        set(gcf, 'color', cfg.bgcol);
    end
    
    function hf = ExamplesFigure(idx)
        % plot example tuning
        
        nsubjects2use = 5;
        
        % ids of extreme cells by metric value
        [~, unselective] = cellfun(@min, SeT(idx));
        [~, selective] = cellfun(@max, SeT(idx));
        [~, unstable] = cellfun(@min, StT(idx));
        [~, stable] = cellfun(@max, StT(idx));
        [~, enhanced] = cellfun(@min, GSS(idx));
        [~, suppressed] = cellfun(@max, GSS(idx));
        
        params.nsubjects2use = nsubjects2use;
        params.idx = idx;
        params.m = m;
        
        hf = figure; row_n = 0;
        
        row_n = row_n+1;
        plotExamples(unselective, Tuning, row_n, cfg, params)
        subplot(6,nsubjects2use,(row_n-1)*nsubjects2use+1)
        ylabel('unselective')
        
        row_n = row_n+1;
        plotExamples(selective, Tuning, row_n, cfg, params)
        subplot(6,nsubjects2use,(row_n-1)*nsubjects2use+1)
        ylabel('selective')
        
        row_n = row_n+1;
        plotExamples(unstable, Tuning, row_n, cfg, params)
        subplot(6,nsubjects2use,(row_n-1)*nsubjects2use+1)
        ylabel('unstable')
        
        row_n = row_n+1;
        plotExamples(stable, Tuning, row_n, cfg, params)
        subplot(6,nsubjects2use,(row_n-1)*nsubjects2use+1)
        ylabel('stable')
        
        row_n = row_n+1;
        plotExamples(enhanced, Tuning, row_n, cfg, params)
        subplot(6,nsubjects2use,(row_n-1)*nsubjects2use+1)
        ylabel('enhanced')
        
        row_n = row_n+1;
        plotExamples(suppressed, Tuning, row_n, cfg, params)
        subplot(6,nsubjects2use,(row_n-1)*nsubjects2use+1)
        ylabel('suppressed')
        
        set(gcf, 'color', cfg.bgcol);
    end

end

%% functions

function plotExamples(cellsArray, metricArray, row_n, cfg,  params)
for i = 1:params.nsubjects2use
    subplot(6,params.nsubjects2use,(row_n-1)*params.nsubjects2use+i)

    thissj = params.idx(i);
    thiscell = cellsArray(i);
    thismat = squeeze(metricArray{thissj}(thiscell,:,:)); % [stims x reps]
    [nStims, nReps] = size(thismat);
    stimlabs = unique(params.m.labels{thissj});

    imagesc(thismat)
    title(['sj#',num2str(thissj),', unit#',num2str(thiscell)])
    xticks(1:nReps); xlabel('repetition')
    yticks(1:nStims); yticklabels(stimlabs)

    colormap(cfg.colormapName)
    axis equal; axis tight
    set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol);
    colorbar
end
end

function plotScatter(metric1, metric2, cfg, color_by)
nsubjects = numel(metric1);

c = defcolor();

metric1 = cell2mat(metric1); % [modes*subjects x 1]
metric2 = cell2mat(metric2); % [modes*subjects x 1]

scatter(metric1,metric2,20,cfg.c(c,:),"filled")

    function c = defcolor()
        c = [];
        for i = 1:nsubjects
            
            switch color_by
                case 'subject'
                    thisc = i*ones(numel(metric1{i}),1);
                otherwise
                    error('coloring metric unknown')
            end
            
            c = [c;thisc];
        end
    end
end

function all_out = loopMetricExtraction(events,labels,method)
nsubjects = numel(events);
all_out = cell(nsubjects,1);
for i = 1:nsubjects
    thisevents = events{i};
    thislabs = labels{i};
    
    % if there are no allowed stimuli here, skip subject
    if isempty(thislabs); continue; end
    
    % call post-processing function
    all_out{i} = extractActivityMetric(thisevents, ...
        method,'cells', ...
        'StimTypes', thislabs);
end
end

