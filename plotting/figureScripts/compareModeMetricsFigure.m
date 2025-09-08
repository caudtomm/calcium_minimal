function hf = compareModeMetricsFigure(v)
    cfg = v.plotConfig;
    
    v.dataFilter.mode_name = 'native_units';
    m = ModeSelector(v).extract; % extract modes
    [GSS, SeT, StT, Tuning] = getSelectedMetrics(m);
    hf = ScatterFigure();
    
    % pick 5 random subjects
    nsubjects = numel(m.data);
    nsubjects2use = 5;
    idx = randi(nsubjects,nsubjects2use,1);
    
    hf = SJwiseCorrelations();
    hf = ExamplesFigure(idx);
    
    % NMF
    v.dataFilter.mode_name = 'nmf';
    v.dataFilter.mode_params.nfactors = 15;
    m = ModeSelector(v).extract; % extract modes
    [GSS, SeT, StT, Tuning] = getSelectedMetrics(m);
    hf = ScatterFigure();
    hf = SJwiseCorrelations();
    hf = ExamplesFigure(idx);
    
    % PCA
    v.dataFilter.mode_name = 'pca';
    v.dataFilter.mode_params.nfactors = 15;
    m = ModeSelector(v).extract; % extract modes
    [GSS, SeT, StT, Tuning] = getSelectedMetrics(m);
    hf = ScatterFigure();
    hf = SJwiseCorrelations();
    hf = ExamplesFigure(idx);

    % dPCA
    v.dataFilter.mode_name = 'dpca';
    v.dataFilter.mode_OI = 'stimulus';
    v.dataFilter.mode_method = 'mode_values';
    m = ModeSelector(v).extract; % extract modes
    [GSS, SeT, StT, Tuning] = getSelectedMetrics(m);
    hf = ScatterFigure();
    hf = SJwiseCorrelations();
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
        metric3 = cell2mat(Tuning);
        metric3 = mean(metric3,[2 3],'omitmissing');
        
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
        plotScatter(StT,SeT,cfg,'subject', metric3)
        box off
        axis square
        xlabel('Stability of Tuning')
        ylabel('Selectivity of Tuning')
        zlabel('Mean activity')
        
        for i = 1:3
            subplot(1,3,i)
            set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol);
        end
        set(gcf, 'color', cfg.bgcol);
    end

    function hf = SJwiseCorrelations()
        hf = figure;

        GSS_StT_corr = zeros(nsubjects,1);
        GSS_SeT_corr = zeros(nsubjects,1);
        SeT_StT_corr = zeros(nsubjects,1);

        for i = 1:nsubjects
            GSS_StT_corr(i) = 1-pdist([GSS{i},StT{i}]', 'correlation');
            GSS_SeT_corr(i) = 1-pdist([GSS{i},SeT{i}]', 'correlation');
            SeT_StT_corr(i) = 1-pdist([SeT{i},StT{i}]', 'correlation');
        end

        y = [GSS_StT_corr,GSS_SeT_corr,SeT_StT_corr];

        boxplot(y)

        xticks(1:3)
        xticklabels({'g. Suppression vs Stability', ...
            'g. Suppression vs Selectivity', ...
            'Selectivity vs Stability'})
        ylabel('correlation')

        box off

        set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol);
        set(gcf, 'color', cfg.bgcol);
        set(gcf, "Position", [0 200 200 500])

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

function plotScatter(metric1, metric2, cfg, color_by, metric3)
nsubjects = numel(metric1);

c = defcolor();

metric1 = cell2mat(metric1); % [modes*subjects x 1]
metric2 = cell2mat(metric2); % [modes*subjects x 1]


if exist("metric3","var") 
    scatter3(metric1,metric2,metric3,20,cfg.c(c,:),"filled")
else
    scatter(metric1,metric2,20,cfg.c(c,:),"filled")
end
% plotHeatmapAndIsoclines(metric1,metric2,20,1,1,1)
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

