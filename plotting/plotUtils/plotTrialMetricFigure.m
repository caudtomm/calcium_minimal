function [hf, out, groups, odor_sets] = plotTrialMetricFigure(v, ps_lim, method, distr_over)  
    arguments
        v ExperimentViewer
        ps_lim = [1 20]
        method = 'population sparseness'
        distr_over = 'cells'
    end
    
    % Knobs
    groups = {'naïve', ...
                'trained', ...
                'trained1', ...
                'trained2', ...
                'trained1 (T-R-S-H-A-ACSF/L)', ...
                'uncoupled'};
    odor_sets = {'all stimuli', ...
                };

    % useful metrics
    ngroups = numel(groups);
    nodor_sets = numel(odor_sets);
    nplots = ngroups*nodor_sets;

    % Initialize output
    hf = gobjects(3,1);
    i_hf = 0;
    out = cell(nplots,1);


    %% Figure 1: comparison of trial metric distributions
    i_hf = i_hf+1;
    plotType = 'boxplot';
    plotGrid;

    %% Figure 2: plot by repetition number
    i_hf = i_hf+1;
    plotType = 'boxplot_repetitions';
    [data,~,n_subplots] = plotGrid;
    set(gcf,'Position',[1 1 300 1000])

    %% Figure 3: imagesc of group and repetition averages
    i_hf = i_hf+1;
    hf(i_hf) = figure; % imagesc plot
    cfg = v.plotConfig;
    
    % Take average along dim 1 of each element in 'data' and concatenate
    avg_data = cellfun(@(x) mean(x,1,'omitmissing'), data, 'UniformOutput', false);
    avg_matrix = vertcat(avg_data{:});
    imagesc(avg_matrix, 'AlphaData', ~isnan(avg_matrix)); % plot with NaNs transparent

    % labels should be stimulus groups by default, and subject groups if there
    % is only one column (=stimulus group)
    labs = repelem(odor_sets,ngroups);
    if nodor_sets == 1
        labs = groups;
    end

    % Cosmetics and labels
    crange = [min(avg_matrix,[],'all','omitmissing'), ...
        max(avg_matrix,[],'all','omitmissing')];
    n_repetitions = size(avg_matrix, 2);
    clim(crange)
    colormap(cfg.colormapName)
    a = colorbar('Color',cfg.axcol);
    a.Label.String = method;
    a.Label.FontSize= gca().FontSize;
    axis equal
    axis tight
    set(gca, 'XTick', 1:n_repetitions, ...
        'YTick', 1:n_subplots, 'YTickLabel', labs);
    xlabel('repetition number')
    set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
    set(gcf, 'color', cfg.bgcol); 
    hold off

    %% Figure 4: comparison of group averages

    i_hf = i_hf+1;
    hf(i_hf) = figure; % single comparison plot
    cfg = v.plotConfig;

    % linearized output distributions (all trials and subjects)
    ytmp = cellfun(@(x) x(:),data,'UniformOutput',false);

    % combine into 1D vectors [labels] , [values]
    labid = [];
    y = [];
    for i = 1:n_subplots
        thisn = length(ytmp{i});
        labid = [labid; repelem(i,thisn,1)];
        y = [y; ytmp{i}];
    end
    
    % labels are inherited from figure 3

    boxplot(y,labid,'Labels',labs,'PlotStyle','compact')

    box off
    ylabel(method)
    set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol);
    set(gcf, 'color', cfg.bgcol); 
    set(gcf,'Position',[1 1 300 1000])
    hold off


    
    %% functions 

    function [out, labs, n_subplots] = plotGrid()

    % define figure size
    ncols = nodor_sets;
    nrows = ngroups;

    hf(i_hf) = figure; % [groups, odor_sets]
    n = 1;
    labs = cell(nplots,1);
    for i_g = 1:ngroups
        % filter data by group
        thisgroup = groups{i_g};
        v.dataFilter.subjectGroup = thisgroup;

        for i_o = 1:nodor_sets
            thisodorset = odor_sets{i_o};

            thisodorset_str = thisodorset;
            if iscell(thisodorset_str); thisodorset_str = strjoin(thisodorset, ', '); end
            msg = ['Plotting group ''',thisgroup,''' for odors: ',thisodorset_str];
            disp(msg)

            % build axes
            subplot(nrows,ncols,n);

            % call intermediate-level plotter
            out{n} = v.plotTrialActivityMetricHead('ps_lim',ps_lim, ...
                            'plotType', plotType, ...
                            'method', method, ...
                            'n_equals', distr_over, ...
                            'trial_sorting', 'stim_id', ...
                            'stim_allowed', thisodorset, ...
                            'reps_touse', [], ...
                            'do_normalize', false);

            % override title and ylabel
            title(thisodorset_str)
            ylabel(thisgroup)

            % if there is no data, delete the subplot
            if isempty(out{n}); axis off; end

            % export label
            labs{n} = [thisgroup,' - ',thisodorset_str];

            % advance axis counter
            n = n+1;
        end
    end

    % set y-axis limits according to globally lowest and highest values
    n_subplots = n-1;

    % get global y-limits
    lims = nan(n_subplots,2); % [min, max]
    for i = 1:n_subplots
        if isempty(out{i}); continue; end

        lims(i,1) = min(out{i},[],'all','omitmissing');
        lims(i,2) = max(out{i},[],'all','omitmissing');
    end
    global_lim = [min(lims(:,1),[],1,'omitmissing'), max(lims(:,2),[],1,'omitmissing')];

    % apply y-limits
    for i = 1:n_subplots
        subplot(nrows,ncols,i)
        if isempty(out{i}); continue; end
        ylim(global_lim)
    end

    set(gcf,'Position',[1 1 2000 1000])

    end
end
