

function [hf, cell_metrics] = compareTraceTypes(v)
    % compare dFoverF and pSpike
    arguments
        v ExperimentViewer
    end 

    % knobs
    nbins = 30; % for deviation from linearity (applied to pSpike)
                % and avg transfer function curve (applied to dF)
    ncellscatterplots = 3;
    sort_bygroup = true;
    
    % init vars
    nfish = numel(v.traces);
    hf = gobjects(9,1);
    cfg = v.plotConfig;
    
    % figure init
    hf(1) = figure;
    ncols = floor(sqrt(nfish))+1;
    nrows = ncols-1;
    
    cell_metrics = cell(nfish, 1);
    dF_all = [];
    pSpike_all = [];
    for i = 1:nfish
        disp(['fish #',num2str(i)])
        
        traces = v.traces{i};
        dF = traces.format(traces.dFoverF_good);
        pSpike = traces.format(traces.pSpike);
    
        [T,N] = size(dF);
        disp([num2str(T),' timepoints, ',num2str(N),' cells.'])
        
        % pSpike vs dF for each cell 
        cells_toplot = randi(N,ncellscatterplots,1);
        polycoef = zeros(N,2);      % linear fit coefficients
        deviations = zeros(N,1);    % 'deviation from linearity' (see Rupprecht et al. 2025)
        ndatapoints = zeros(N,1);
        normr = zeros(N,1); 
        curves = zeros(nbins,N);
        for i_cell = 1:N
            thisdF = dF(:,i_cell);
            thispSpike = pSpike(:,i_cell);
    
            % clean up nan values
            to_rem = isnan(thisdF) | isnan(thispSpike);
            thisdF = thisdF(~to_rem);
            thispSpike = thispSpike(~to_rem);
            ndatapoints(i_cell) = length(thisdF);
            dF_all = [dF_all; thisdF];
            pSpike_all = [pSpike_all; thispSpike];
    
            % linear fit
            [p,S] = polyfit(thisdF, thispSpike, 1);
            polycoef(i_cell,:) = p;
    
            % (debug) plot single cell scatter
            % if ismember(i_cell,cells_toplot)
            %     figure; scatter(thisdF,thispSpike,4,'k','filled')
            %     axis square tight;
            %     x = xlim; y = polyval(p,x);
            %     hold on; line(x,y,'Color','r','LineWidth',1.5)
            %     xlabel('dFoverF'); ylabel('inferred SR')
            %     title(['cell #',num2str(i_cell)])
            % end
            
            % deviation from linearity
            thispSpike_linear = arrayfun(@(x) p(1)*x+p(2),thisdF);
            bin_edges = linspace(min(thispSpike),max(thispSpike),nbins+1);
            deviations_bins = zeros(nbins,1);
            for i_bin = 1:nbins
                % get average values in the bin
                idx = thispSpike>=bin_edges(i_bin) & thispSpike<=bin_edges(i_bin+1);
                binpSpike = mean(thispSpike(idx));
                binpSpike_linear = mean(thispSpike_linear(idx));
                
                % get deviation
                deviations_bins(i_bin) = (binpSpike - binpSpike_linear)^2 / (binpSpike_linear^2);
            end
            deviations(i_cell) = sqrt(mean(deviations_bins,'omitmissing'));
    
            % norm of fit residuals
            normr(i_cell) = S.normr;
    
            % avg transfer function curve
            bin_edges = linspace(min(dF,[],'all'),max(dF,[],'all'),nbins+1); % bins over all cells
            for i_bin = 1:nbins
                % get average value in the bin
                idx = thisdF>=bin_edges(i_bin) & thisdF<=bin_edges(i_bin+1);
                curves(i_bin,i_cell) = mean(thispSpike(idx));
            end
        end
    
        % store single-cell metrics to table
        slopes = polycoef(:,1);
        intercepts = polycoef(:,2);
        noise_level_PR = traces.dFnoise(traces.goodNeuron_IDs);
        pxVariance_overtime = traces.format(traces.Fnoise(:,traces.goodNeuron_IDs,:));
        pxVariance = mean(pxVariance_overtime,'omitmissing');
        pxVariance = pxVariance(:);
        cell_metrics{i} = table( ...
            ndatapoints, ...
            deviations, ...
            slopes, ...
            intercepts, ...
            normr, ...
            noise_level_PR, ...
            pxVariance);
    
        % transfer function plot
        figure(hf(1));
        subplot(nrows,ncols,i)
        x = repmat(bin_edges(2:end),N,1)'; % [nbins,N]
        plot(x,curves,'k')
        axis square tight; box off
        xlabel('dFoverF'); ylabel('inferred SR');
        title(['fish #',num2str(i)])
    end
    
    % sort order of fish by group
    [sorted_groups,idxbygroup] = sort(v.subjectTab.group);
    cell_metrics_sorted = cell_metrics(idxbygroup);
    grouplabels = sorted_groups;
    
    % cell metrics boxplots
    labels = cell_metrics{1}.Properties.VariableNames;
    for i = 1:numel(labels)
        thislabel = labels{i};
        values = [];
        for i_fish = 1:nfish
            thisvalues = cell_metrics_sorted{i_fish}.(thislabel);
            N = numel(thisvalues);
            values = [values; thisvalues repelem(i_fish,N,1)];
        end
        hf(1+i) = figure;
        boxplot(values(:,1),values(:,2),'PlotStyle','compact','Symbol','');
        xticks(1:nfish); xticklabels(grouplabels);xtickangle(90)
        title(thislabel)
        box off
    end
    
    % heatmap of inferred SR vs dFoverF over all data
    hf(9) = figure;
    subplot(131); h(1)=histogram(dF_all,50);
    xlabel('dFoverF','Color',cfg.textcol);
    ylabel('histogram','Color',cfg.textcol);
    axis square; box off
    subplot(132); h(2)=histogram(pSpike_all,50);
    xlabel('inferred SR','Color',cfg.textcol);
    axis square; box off
    subplot(133); plotHeatmapAndIsoclines(dF_all, pSpike_all, 50, 1,0,1)
    axis square
    b = colorbar; b.Label.String = 'Log Density';
    xlabel('dFoverF','Color',cfg.textcol);
    ylabel('inferred SR','Color',cfg.textcol);
    title('')
    set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol);
    set(gcf, 'color', cfg.bgcol); 
    set(gcf,'Position',[100 100 1000 300])


end