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

    % spearman instead iof pearson corr, because we cannot assume normal distributions
    hf = SJwiseCorrelations('spearman'); 
    hf = ExamplesFigure(idx);
    
    % % NMF
    % v.dataFilter.mode_name = 'nmf';
    % v.dataFilter.mode_params.nfactors = 15;
    % m = ModeSelector(v).extract; % extract modes
    % [GSS, SeT, StT, Tuning] = getSelectedMetrics(m);
    % hf = ScatterFigure();
    % hf = SJwiseCorrelations();
    % hf = ExamplesFigure(idx);
    % 
    % % PCA
    % v.dataFilter.mode_name = 'pca';
    % v.dataFilter.mode_params.nfactors = 15;
    % m = ModeSelector(v).extract; % extract modes
    % [GSS, SeT, StT, Tuning] = getSelectedMetrics(m);
    % hf = ScatterFigure();
    % hf = SJwiseCorrelations();
    % hf = ExamplesFigure(idx);

    % % dPCA
    % v.dataFilter.mode_name = 'dpca';
    % v.dataFilter.mode_OI = 'stimulus';
    % v.dataFilter.mode_method = 'mode_values';
    % m = ModeSelector(v).extract; % extract modes
    % [GSS, SeT, StT, Tuning] = getSelectedMetrics(m);
    % hf = ScatterFigure();
    % hf = SJwiseCorrelations();
    % hf = ExamplesFigure(idx);

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
        metric1 = StT; metric2 = GSS;
        plotScatter(metric1,metric2,cfg,'subject')
        plotHeatmapAndIsoclines(cell2mat(metric1),cell2mat(metric2),20,1,0,0);
        crange = clim;
        plotHeatmapAndIsoclines(cell2mat(metric1),cell2mat(metric2),10,0,1,0);
        clim(crange);
        xlim([-1 1]);
        ylim([-3 3])
        colormap(cfg.colormapName)
        box off
        axis square
        xlabel('Stability of Tuning')
        ylabel('g. Suppression')
        
        subplot(132)
        metric1 = SeT; metric2 = GSS;
        plotScatter(metric1,metric2,cfg,'subject')
        plotHeatmapAndIsoclines(cell2mat(metric1),cell2mat(metric2),20,1,0,0);
        crange = clim;
        plotHeatmapAndIsoclines(cell2mat(metric1),cell2mat(metric2),10,0,1,0);
        xlim([0 1]);
        ylim([-3 3])
        clim(crange);
        colormap(cfg.colormapName)
        box off
        axis square
        xlabel('Selectivity of Tuning')
        ylabel('g. Suppression')
        
        subplot(133)
        % plotScatter(StT,SeT,cfg,'subject', metric3)
        metric1 = StT; metric2 = SeT;
        plotScatter(metric1,metric2,cfg,'subject')
        plotHeatmapAndIsoclines(cell2mat(metric1),cell2mat(metric2),20,1,0,0);
        crange = clim;
        plotHeatmapAndIsoclines(cell2mat(metric1),cell2mat(metric2),10,0,1,0);
        xlim([-1 1]);
        ylim([0 1])
        clim(crange);
        colormap(cfg.colormapName)
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

    function hf = SJwiseCorrelations(method)
        hf = figure;

        GSS_StT_corr = zeros(nsubjects,1);
        GSS_SeT_corr = zeros(nsubjects,1);
        SeT_StT_corr = zeros(nsubjects,1);

        for i = 1:nsubjects
            GSS_StT_corr(i) = 1-pdist([GSS{i},StT{i}]', method);
            GSS_SeT_corr(i) = 1-pdist([GSS{i},SeT{i}]', method);
            SeT_StT_corr(i) = 1-pdist([SeT{i},StT{i}]', method);
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

        % stats over full dataset
        [rho,pval] = corr([cell2mat(GSS),cell2mat(StT)],'Type','Spearman','Rows','pairwise');
        disp(['GSS vs StT: corr=',num2str(rho(2)),', p=',num2str(pval(2))])
        [rho,pval] = corr([cell2mat(GSS),cell2mat(SeT)],'Type','Spearman','Rows','pairwise');
        disp(['GSS vs SeT: corr=',num2str(rho(2)),', p=',num2str(pval(2))])
        [rho,pval] = corr([cell2mat(SeT),cell2mat(StT)],'Type','Spearman','Rows','pairwise');
        disp(['SeT vs StT: corr=',num2str(rho(2)),', p=',num2str(pval(2))])

        % linear mixed-effects model
        [~,r_rm,pVal] = doLMEfit(GSS,StT);
        fprintf('Repeated Measures Correlation (%s-%s): r = %.3f, p = %.4f\n','GSS','StT', r_rm, pVal);
        [~,r_rm,pVal] = doLMEfit(GSS,SeT);
        fprintf('Repeated Measures Correlation (%s-%s): r = %.3f, p = %.4f\n','GSS','SeT', r_rm, pVal);
        [~,r_rm,pVal] = doLMEfit(SeT,StT);
        fprintf('Repeated Measures Correlation (%s-%s): r = %.3f, p = %.4f\n','SeT','StT', r_rm, pVal);

    end

    function [lme,r_rm,pVal] = doLMEfit(metric1,metric2)
        % metric1 = predictor (x), metric2 = response (y), both {nsubjects} cells
        result = statsUtils.lmeRegress(metric2, metric1);
        lme  = result.lme;
        r_rm = result.r_rm;
        pVal = result.p;
    end
    
    function hf = ExamplesFigure(idx)
        % plot example tuning
        
        nsubjects2use = 5;
        
        % ids of extreme cells by metric value
        [vals.unselective, unselective] = cellfun(@min, SeT(idx));
        [vals.selective, selective] = cellfun(@max, SeT(idx));
        [vals.unstable, unstable] = cellfun(@min, StT(idx));
        [vals.stable, stable] = cellfun(@max, StT(idx));
        [vals.enhanced, enhanced] = cellfun(@min, GSS(idx));
        [vals.suppressed, suppressed] = cellfun(@max, GSS(idx));
        
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

        vals
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

