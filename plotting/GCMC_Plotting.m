classdef GCMC_Plotting
methods (Static)
    %% no data required
    
    function hf = coverageFigure(cfg)
        arguments
            cfg PlotConfig = PlotConfig()
        end

        manifold_size_range = [400 1000]; % from this study, after taking out nans
        Bo_manifold_size = 113; % 5s * 7.5 Hz * 3 trials/odor

        n_points_range = [1:10:151];

        minCvg = GCMC_Analysis.estimateManifoldCoverage(manifold_size_range(2), n_points_range);
        maxCvg = GCMC_Analysis.estimateManifoldCoverage(manifold_size_range(1), n_points_range);
        BoCvg = GCMC_Analysis.estimateManifoldCoverage(Bo_manifold_size, n_points_range);

        estCvg = minCvg; % arbitrary

        hf = figure;
        imagesc(estCvg.X,estCvg.Y,maxCvg.expectedCoverage-minCvg.expectedCoverage)
        hold on

        scatter(minCvg.full_coverage,estCvg.Y,50,'g','filled')
        scatter(maxCvg.full_coverage,estCvg.Y,50,'r','filled')
        scatter(BoCvg.full_coverage,estCvg.Y,50,'k','filled')
        
        b(1) = plot(minCvg.full_coverage,estCvg.Y,'g-','LineWidth',1.5);
        b(2) = plot(maxCvg.full_coverage,estCvg.Y,'r-','LineWidth',1.5);
        b(3) = plot(BoCvg.full_coverage,estCvg.Y,'k-','LineWidth',1.5);

        axis equal tight
        xlabel(estCvg.X_label); ylabel(estCvg.Y_label)
        title('expected point cloud coverage (max - min)')
        legend(b,{'min','max','Hu et al. 2024'},'BackgroundAlpha',0,'Box','off','TextColor',cfg.bgcol)
        
        colormap(cfg.colormapName)
        set(gca, 'color', cfg.bgcol, 'XColor', cfg.axcol, 'YColor', cfg.axcol, 'ZColor', cfg.axcol);
        set(gcf, 'color', cfg.bgcol);

    end

    function params = manifold_size2params_Map(size_range, min_N)
        % knobs

        % point num sampling
        points_stepsz = 10; % step size
        points_lim = [1 101]; 
        % arbitrarily add 40 repetitions, to ensure full coverage. empirically,
        % metric error curves saturate (go sublinear in loglog space) fairly late compared to
        % the full coverage frontier, so we want to be generous. oversampling is not an issue
        % here, except for computational demands. offset will scale compute requirements linearly.
        frontier_offset = 40;

        % usefuls
        n_sizes = numel(size_range);
        max_points = floor(min_N/2)-5; % the 5 is an arbitrary offset
        points_range = points_lim(1):points_stepsz:points_lim(2);

        % for each size, find optima parameters
        params = zeros(n_sizes,3);
        for i = 1:n_sizes
            this_size = size_range(i);
            estCov = GCMC_Analysis.estimateManifoldCoverage(this_size,points_range);
            
            frontier = estCov.full_coverage;

            size_tolerance = floor(.3 * this_size);
            best_npoints = min([max_points, this_size-size_tolerance]);
            [~,idx] = min(abs(points_range-best_npoints));
            best_repnum = frontier_offset + frontier(idx);

            params(i,:) = [this_size, best_npoints, best_repnum];
        end
    end

    function results = estimateManifoldCoverage(manifold_size,n_points)
        arguments
            manifold_size double {mustBePositive}
            n_points double {mustBePositive}
        end

        % analytical estimate of manifold coverage based on number of points sampled
        max_repetitions = 300;
        n_random_samples = 100;
        bin_size = 10; % 10 repetitions increments
        nBins = floor(max_repetitions / bin_size);

        bins = bin_size * (1:nBins);
        n_n_points = numel(n_points);
        avg_curve = zeros(n_n_points,nBins);
        for i_points = 1:n_n_points
            this_n_points = n_points(i_points);
            fprintf('Estimating coverage for samples of %s points\n',num2str(this_n_points));
            for i_bin = 1:nBins
                n_reps = bins(i_bin);
                coverage_vals = zeros(1,n_random_samples); % 100 random samples
                for i_sample = 1:n_random_samples
                    % simulate point sampling
                    points = randi(manifold_size, this_n_points, n_reps);
                    coverage_vals(i_sample) = numel(unique(points)) / manifold_size;
                end
                avg_curve(i_points,i_bin) = mean(coverage_vals);

            end
        end

        % find first x-index (num of repetitions) that reach full
        % coverage for each point cloud size
        full_coverage = nBins*ones(n_n_points,1);
        for i = 1:n_n_points
            val = find(avg_curve(i,:)>=0.99,1,'first');
            if ~isempty(val); full_coverage(i) = val; end
        end

        results.expectedCoverage = avg_curve;
        results.X = bins;
        results.X_label = 'repetition number';
        results.Y = n_points;
        results.Y_label = 'size of point cloud';
        results.manifold_size = manifold_size;
        results.n_points = n_points;
        results.n_random_samples = n_random_samples;
        results.full_coverage = bins(full_coverage);

    end
    
    %% yes data required

    function [hf, all_stats] = plotBoxplotsByGroup(group_data, cfg, filter, metrics)
        arguments
            group_data struct
            cfg PlotConfig = PlotConfig()
            filter GCMCResultsFilter = GCMCResultsFilter('shuffle', false)
            metrics cell = {}
        end

        nGroups = numel(group_data);
        group_names = {group_data.group_name};
        if isempty(metrics)
            metrics = group_data(1).data.Properties.VariableNames(5:end-3);
        end
        nMetrics = numel(metrics);

        hf = gobjects(nMetrics,1);
        all_stats = struct('metric', {}, 'mw', {});

        for i_metric = 1:nMetrics
            hf(i_metric) = figure;
            hold on;

            disp(['Metric: ', metrics{i_metric}])

            % Collect data for each group
            datacells = cell(1,nGroups);
            for i_group = 1:nGroups
                filtered_table = filter.filterTable(group_data(i_group).data);
                this_data = filtered_table{:, metrics{i_metric}};
                datacells{i_group} = this_data(:);
            end

            % Pairwise Mann-Whitney U-tests with FDR correction
            mw_results = statsUtils.pairwiseMannWhitney(datacells, group_names);
            all_stats(i_metric).metric = metrics{i_metric};
            all_stats(i_metric).mw = mw_results;

            for k = 1:numel(mw_results)
                disp(sprintf('  %s vs %s: p_raw=%.4g, p_adj=%.4g %s', ...
                    mw_results(k).group1, mw_results(k).group2, ...
                    mw_results(k).p_raw, mw_results(k).p_adj, ...
                    statsUtils.pvalToStars(mw_results(k).p_adj)));
            end

            % Create boxplot
            RF_mkBoxPlot3(datacells,[],[],.1,.5,.5,2,[]);
            xticks(1:nGroups); xticklabels(group_names)

            title(['Metric: ', metrics{i_metric}]);
            ylabel(metrics{i_metric});

            % Customize plot appearance
            box off;
            axis tight
            xlim([.5 max(xticks)+.5])

            % Add significance brackets
            yl = ylim;
            y_range = yl(2) - yl(1);
            bracket_y = yl(2) + 0.05 * y_range;
            bracket_step = 0.08 * y_range;
            for k = 1:numel(mw_results)
                if mw_results(k).p_adj < 0.05
                    addSignificanceAnnotation(gca, ...
                        [mw_results(k).i, mw_results(k).j], ...
                        bracket_y, mw_results(k).p_adj, ...
                        'style', 'bracket', 'color', cfg.axcol);
                    bracket_y = bracket_y + bracket_step;
                end
            end
            ylim([yl(1), bracket_y + 0.02 * y_range]);

            set(gca, 'color', cfg.bgcol, 'XColor', cfg.axcol, 'YColor', cfg.axcol, 'ZColor', cfg.axcol);
            set(gcf, 'color', cfg.bgcol);
            set(gcf, 'Position', [100, 100, 200, 500]);
            hold off;
        end

    end

    function plotBoxplotsForEachMetric(avg_results, cfg, filter, metrics)
        arguments
            avg_results table % from a single subj
            cfg PlotConfig = PlotConfig()
            filter GCMCResultsFilter = GCMCResultsFilter()
            metrics cell = {}
        end

        if isempty(metrics)
            metrics = avg_results.Properties.VariableNames(5:end-2);
        end

        % Apply filter (excluding shuffle filter, which is handled separately below)
        filter_no_shuffle = filter;
        filter_no_shuffle.shuffle = [];  % don't filter by shuffle here
        filtered_results = filter_no_shuffle.filterTable(avg_results);

        % separate data and shuffle conditions
        y_data = filtered_results{filtered_results.shuffle==false, metrics};
        y_shuffle = filtered_results{filtered_results.shuffle==true, metrics};

        % create boxplots for each metric
        nMetrics = numel(metrics);
        for i = 1:nMetrics
            figure;
            hold on;

            thisdata = y_data(:,i);
            thisshuffle = y_shuffle(:,i);

            % Mann-Whitney U-test (non-parametric, unpaired)
            p_mw = NaN;
            if numel(thisdata) >= 2 && numel(thisshuffle) >= 2
                p_mw = ranksum(thisdata, thisshuffle);
            end

            % Paired t-test (one-tailed, using mean difference to pick tail)
            mean_data = mean(thisdata);
            mean_shuffle = mean(thisshuffle);
            if mean_data > mean_shuffle; tail = 'right'; else; tail = 'left'; end
            try
                [~, p_tt] = ttest(thisdata, thisshuffle, 'Tail', tail);
            catch
                p_tt = NaN;
            end

            disp(sprintf('Metric: %s | Mann-Whitney p=%.4g | paired t-test p=%.4g', ...
                metrics{i}, p_mw, p_tt));

            % Box plot
            datacells = {thisdata, thisshuffle};
            RF_mkBoxPlot3(datacells,[],[],.1,.5,.5,2,[]);
            xticks([1 2]); xticklabels({'Data', 'Shuffle'})

            % Scatter overlay
            scatter(repelem(1, size(y_data, 1)), y_data(:, i), 'r', 'filled', 'jitter', 'on', 'jitterAmount', 0.15);
            scatter(repelem(2, size(y_shuffle, 1)), y_shuffle(:, i), 'b', 'filled', 'jitter', 'on', 'jitterAmount', 0.15);

            % Significance bracket (use Mann-Whitney p-value)
            yl = ylim;
            bracket_y = yl(2) + 0.05 * (yl(2) - yl(1));
            addSignificanceAnnotation(gca, [1 2], bracket_y, p_mw, ...
                'style', 'bracket', 'color', cfg.axcol);

            title(metrics{i});
            ylabel(metrics{i});
            box off
            set(gcf,'Position',[100 100 200 400])
            set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
            set(gcf, 'color', cfg.bgcol);
            hold off;
        end
    end

    function plotMetricStability(results, cfg, filter)
        arguments
            results table
            cfg PlotConfig = PlotConfig()
            filter GCMCResultsFilter = GCMCResultsFilter('shuffle', false)
        end

        results = filter.filterTable(results);

        bin_size = 10; % number of measurements per bin
        nRandomSamples = 10; % number of random sampling for each bin
        avg_results = genResultAvgTable(results, []); % stimulus names not needed here
        metrics = avg_results.Properties.VariableNames(5:11); % exclude grouping variables and stimulus names
        N = height(avg_results);
        nMeasurements = avg_results.GroupCount(1); % assume all groups have the same number of measurements
        nBins = floor(nMeasurements / bin_size);

        % get averaged results for incremental repetition bins (each consecutive bin adds bin_size more measurement reps)
        nMetrics = numel(metrics);
        binned_vals = zeros(N, nMetrics, nBins, nRandomSamples);
        for i_bin = 1:nBins
            for i_sample = 1:nRandomSamples
                % randomly sample bin_size * i_bin measurements from each group
                idx = false(nMeasurements,1);
                idx(randperm(nMeasurements, bin_size * i_bin)) = true;

                res = genResultAvgTable(results, idx);

                binned_vals(:,:,i_bin,i_sample) = table2array(res(:,metrics)); % stimulus names not needed here
            end
        end

        % get coefficient of variation across random samples for each bin
        cv_vals = zeros(N, nMetrics, nBins);
        for i_bin = 1:nBins
            for i_metric = 1:nMetrics
                data = squeeze(binned_vals(:,i_metric,i_bin,:)); % N x nRandomSamples
                cv_vals(:,i_metric,i_bin) = abs(std(data,0,2) ./ mean(data,2));
            end
        end

        % plot
        figure;
        hold on;
        mean_handles = gobjects(1, nMetrics); % Preallocate for legend handles
        for i_metric = 1:nMetrics
            % Calculate mean and standard deviation across stimuli
            mean_curve = mean(squeeze(cv_vals(:, i_metric, :)), 1);
            std_curve = std(squeeze(cv_vals(:, i_metric, :)), 0, 1);

            % Plot mean curve
            mean_handles(i_metric) = plot(1:nBins, mean_curve, '-o', 'DisplayName', metrics{i_metric});

            % Add shaded area for standard deviation
            fill([1:nBins, fliplr(1:nBins)], ...
                    [mean_curve + std_curve, fliplr(mean_curve - std_curve)], ...
                    cfg.c(i_metric, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        end
        hold off;

        xticks(1:nBins);
        xticklabels(bin_size * (1:nBins));
        xlabel('Number of Measurements Used');
        ylabel('Coefficient of Variation over random samples');
        title('Stability of Metrics');
        legend(mean_handles, metrics, 'Location', 'best','Box','off'); % Only label mean curves
        box off;
        axis tight
        set(gcf, 'Position', [100 100 600 400]);
        set(gca, 'color', cfg.bgcol, 'XColor', cfg.axcol, 'YColor', cfg.axcol, 'ZColor', cfg.axcol);
        set(gcf, 'color', cfg.bgcol);

    end

    function [hf, stats] = plotMetricsByRepLine(group_data, cfg, filter, metric_name, opts)
        % Plot metric values across repetitions as line plot with error ribbons
        %
        % Uses statsUtils for statistics and plotLineNShade for visualization.
        %
        % Inputs:
        %   group_data - struct array with .group_name and .data (table)
        %   cfg - PlotConfig object
        %   filter - GCMCResultsFilter
        %   metric_name - name of metric column (e.g., 'Fun_capacity')
        %   opts.average_by_subject - if true, average rows per subject before stats
        %
        % Outputs:
        %   hf - figure handle
        %   stats - struct with .friedman and .between_group results

        arguments
            group_data struct
            cfg PlotConfig = PlotConfig()
            filter GCMCResultsFilter = GCMCResultsFilter('shuffle', false)
            metric_name char = 'Fun_capacity'
            opts.average_by_subject logical = false
        end

        nGroups = numel(group_data);
        group_names = {group_data.group_name};

        % === Extract data by group and rep ===
        [data_by_group_rep, reps] = extractDataByGroupRep(group_data, filter, metric_name, opts.average_by_subject);
        nReps = numel(reps);

        if nReps == 0
            warning('No repetition data found.');
            hf = figure; stats = struct();
            return;
        end

        % === Compute summary statistics ===
        [means, sems] = computeGroupRepStats(data_by_group_rep);

        % === Statistical tests ===
        stats = struct();

        % Friedman test per group
        stats.friedman = computeFriedmanByGroup(group_data, filter, metric_name, reps);

        % Between-group Mann-Whitney at each rep (with FDR)
        stats.between_group = computeBetweenGroupStats(data_by_group_rep, group_names, reps);

        % === Plotting ===
        hf = figure;
        hold on;
        line_handles = gobjects(nGroups, 1);

        for i_group = 1:nGroups
            color = cfg.c(i_group, :);

            % Use plotLineNShade for ribbon
            mu = means(i_group, :);
            sem = sems(i_group, :);
            valid = ~isnan(mu);

            if sum(valid) > 1
                lcfg.lineStyle = '-';
                lcfg.lineWidth = .5;
                plotLineNShade(reps(valid), mu(valid)', sem(valid)', color, lcfg);
            end

            % Build legend with Friedman result
            legend_str = group_names{i_group};
            if i_group <= numel(stats.friedman) && ~isnan(stats.friedman(i_group).p)
                legend_str = sprintf('%s (F: %s)', legend_str, ...
                    statsUtils.pvalToStars(stats.friedman(i_group).p));
            end

            line_handles(i_group) = plot(reps, mu, 'o', ...
                'Color', color, 'MarkerFaceColor', color, ...
                'MarkerSize', 1, 'DisplayName', legend_str);
        end

        % Add significance stars for between-group comparisons
        y_max = max(means(:) + sems(:), [], 'omitnan');
        y_min = min(means(:) - sems(:), [], 'omitnan');
        star_offset = 0.05 * (y_max - y_min);

        for i = 1:numel(stats.between_group)
            if stats.between_group(i).p_fdr < 0.05
                i_rep = stats.between_group(i).i_rep;
                y_pos = y_max + star_offset * (1 + mod(i-1, 3));  % stagger if multiple
                addSignificanceAnnotation(gca, reps(i_rep), y_pos, ...
                    stats.between_group(i).p_fdr, 'color', cfg.axcol);
            end
        end

        hold off;

        % Formatting
        xlabel('Repetition');
        ylabel(strrep(metric_name, 'Fun_', ''));
        legend(line_handles, 'Location', 'best', 'Box', 'off');
        box off;
        xlim([min(reps)-0.5, max(reps)+0.5]);
        xticks(reps);

        set(gca, 'color', cfg.bgcol, 'XColor', cfg.axcol, 'YColor', cfg.axcol);
        set(gcf, 'color', cfg.bgcol);
        set(gcf, 'Position', [100, 100, 400, 350]);
    end

    %% Helper methods

    function [group_data, all_results] = loadAndPrepareData(v, folder_tag, exp_name)
        % Load GCMC results and prepare group data with merged trained groups
        %
        % Inputs:
        %   v          - ExperimentViewer object
        %   folder_tag - subfolder name (e.g., 'odors', 'trials')
        %   exp_name   - experiment name (subfolder under manifold_data/folder_tag)
        %
        % Outputs:
        %   group_data   - struct array with .group_name and .data
        %   all_results  - raw results table

        arguments
            v
            folder_tag char
            exp_name char = ''
        end

        if isempty(exp_name)
            % Try to derive from experiment record
            if isfield(v.experiment.record, 'tabpath')
                [~, exp_name] = fileparts(v.experiment.record.tabpath);
            else
                error('GCMC_Plotting:loadAndPrepareData', ...
                    'exp_name must be specified when experiment record has no tabpath.');
            end
        end

        indir = fullfiletol('manifold_data', folder_tag);
        all_results = GCMC_Analysis(v).extractResultsFromMultipleSubjects(fullfiletol(indir, exp_name));
        group_data_raw = GCMC_Analysis(v).clusterByGroup(all_results);
        group_data = GCMC_Plotting.mergeTrainedGroups(group_data_raw);
    end

    function new_group_data = mergeTrainedGroups(group_data)
        % Merge trained subgroups into single 'trained' group
        new_group_data = struct('group_name', {}, 'data', {});

        naive_idx = find(strcmp({group_data.group_name}, 'naïve'));
        uncoupled_idx = find(strcmp({group_data.group_name}, 'uncoupled'));
        trained_idx = find(contains({group_data.group_name}, 'trained'));

        idx = 0;
        if ~isempty(naive_idx)
            idx = idx + 1;
            new_group_data(idx) = group_data(naive_idx);
        end
        if ~isempty(trained_idx)
            idx = idx + 1;
            new_group_data(idx).group_name = 'trained';
            new_group_data(idx).data = group_data(trained_idx(1)).data;
            for i = 2:numel(trained_idx)
                new_group_data(idx).data = [new_group_data(idx).data; group_data(trained_idx(i)).data];
            end
        end
        if ~isempty(uncoupled_idx)
            idx = idx + 1;
            new_group_data(idx) = group_data(uncoupled_idx);
        end
    end

    function plotAndSaveBoxplots(group_data, cfg, filter, metrics, labels, yranges, savedir, saveType)
        % Plot and save boxplots for multiple metrics
        %
        % Inputs:
        %   group_data - struct array with .group_name and .data
        %   cfg        - PlotConfig
        %   filter     - GCMCResultsFilter
        %   metrics    - cell array of metric column names
        %   labels     - cell array of display labels
        %   yranges    - cell array of [ymin ymax] ranges
        %   savedir    - directory to save figures
        %   saveType   - 'vector' or 'raster'

        hf = GCMC_Plotting.plotBoxplotsByGroup(group_data, cfg, filter, metrics);

        nFigs = numel(hf);
        for i = 1:nFigs
            if i <= numel(labels) && i <= numel(yranges)
                figure(hf(i));
                ylim(yranges{i});
                title('');
                ylabel(labels{i});
                cfg.figSize = 'tiny';
                cfg.aspRatioType = 'tall';
                cfg.lineWidth = 0.5;
                cfg.setLines = true;
                cfg.setFigure;
                cfg.savePath = savedir;
                cfg.saveFigure(gcf, [labels{i}, ' boxplot'], saveType);
            end
        end
        close all;
    end

    function plotAndSaveRepLines(group_data, cfg, filter, metrics, labels, savedir, saveType, opts)
        % Plot and save line plots across repetitions for multiple metrics
        %
        % Inputs:
        %   group_data - struct array with .group_name and .data
        %   cfg        - PlotConfig
        %   filter     - GCMCResultsFilter
        %   metrics    - cell array of metric column names
        %   labels     - cell array of display labels
        %   savedir    - directory to save figures
        %   saveType   - 'vector' or 'raster'
        %   opts.average_by_subject - if true, average rows per subject before stats

        arguments
            group_data struct
            cfg PlotConfig
            filter GCMCResultsFilter
            metrics cell
            labels cell
            savedir
            saveType
            opts.average_by_subject logical = false
        end

        for i = 1:numel(metrics)
            [hf, stats] = GCMC_Plotting.plotMetricsByRepLine(group_data, cfg, filter, metrics{i}, ...
                'average_by_subject', opts.average_by_subject);

            % Display stats summary
            disp(['=== Stats for ', labels{i}, ' ===']);
            for j = 1:numel(stats.friedman)
                disp(sprintf('  Friedman %s: p=%.4g', stats.friedman(j).group, stats.friedman(j).p));
            end
            sig_between = find([stats.between_group.p_fdr] < 0.05);
            for j = sig_between
                disp(sprintf('  Between-group rep%d %s vs %s: p_fdr=%.4g', ...
                    stats.between_group(j).rep, stats.between_group(j).group1, ...
                    stats.between_group(j).group2, stats.between_group(j).p_fdr));
            end

            ylabel(labels{i});
            title('');
            cfg.figSize = 'small';
            cfg.aspRatioType = 'square';
            cfg.setFigure;
            cfg.savePath = savedir;
            cfg.saveFigure(gcf, [labels{i}, ' by rep'], saveType);
            close(hf);
        end
    end

    %% Sliding window methods

    function sw = loadSlidingWindowData(v, folder_tag, exp_name)
        % Load GCMC results across time windows
        %
        % Directory structure:
        %   manifold_data/<folder_tag>/<exp_name>/<time_window>/manifolds_subj*/...
        %   where time_window is e.g. '0.5-4.5s', '-4-0s', '-6--2s'
        %
        % Inputs:
        %   v          - ExperimentViewer object
        %   folder_tag - subfolder name (e.g., 'odors_slide_windows')
        %   exp_name   - experiment name (subfolder under folder_tag)
        %
        % Output:
        %   sw struct with fields:
        %     .windows   - cell{1,nW} of sorted window label strings
        %     .t_centers - [1 x nW] double, window midpoints (seconds)
        %     .t_edges   - [nW x 2] double, [t0, t1] per window
        %     .groups    - cell{1,nG} of group names
        %     .data      - {nG x nW} cell array of tables

        arguments
            v
            folder_tag char
            exp_name char = ''
        end

        if isempty(exp_name)
            if isfield(v.experiment.record, 'tabpath')
                [~, exp_name] = fileparts(v.experiment.record.tabpath);
            else
                error('GCMC_Plotting:loadSlidingWindowData', ...
                    'exp_name must be specified when experiment record has no tabpath.');
            end
        end

        basedir = fullfiletol('manifold_data', folder_tag, exp_name);

        % Discover subdirectories
        entries = dir(basedir);
        entries = entries([entries.isdir] & ~ismember({entries.name}, {'.', '..'}));

        % Filter to valid time-window names and parse edges
        all_names = {entries.name};
        [t_edges_raw, valid] = parseSlidingWindowEdges(all_names);
        window_labels = all_names(valid);
        t_edges_raw = t_edges_raw(valid, :);

        if isempty(window_labels)
            error('GCMC_Plotting:loadSlidingWindowData', ...
                'No valid time window directories found in %s', basedir);
        end

        % Sort by window start time
        [t_edges, sort_idx] = sortrows(t_edges_raw, 1);
        window_labels = window_labels(sort_idx);
        t_centers = mean(t_edges, 2)';
        nWindows = numel(window_labels);

        % Load GCMC results for each window
        gcmc = GCMC_Analysis(v);
        gd_per_window = cell(nWindows, 1);
        for i_win = 1:nWindows
            winpath = fullfiletol(basedir, window_labels{i_win});
            fprintf('Loading window %d/%d: %s\n', i_win, nWindows, window_labels{i_win});
            results = gcmc.extractResultsFromMultipleSubjects(winpath);
            gd_raw = gcmc.clusterByGroup(results);
            gd_per_window{i_win} = GCMC_Plotting.mergeTrainedGroups(gd_raw);
        end

        % Determine consistent group ordering from first window
        group_names = {gd_per_window{1}.group_name};
        nGroups = numel(group_names);

        % Organize into {nGroups x nWindows} cell array of tables
        data = cell(nGroups, nWindows);
        for i_win = 1:nWindows
            gd = gd_per_window{i_win};
            for j = 1:numel(gd)
                grp_idx = find(strcmp(group_names, gd(j).group_name), 1);
                if ~isempty(grp_idx)
                    data{grp_idx, i_win} = gd(j).data;
                end
            end
        end

        sw.windows = window_labels;
        sw.t_centers = t_centers;
        sw.t_edges = t_edges;
        sw.groups = group_names;
        sw.data = data;
    end

    function [hf, p_friedman] = plotSlidingWindowMetrics(sw, cfg, filter, metrics, opts)
        % Plot GCMC metrics across sliding time windows
        %
        % Auto-selects plot mode:
        %   nReps > 1 → imagesc heatmap [reps x windows] per (group, metric)
        %   nReps ≤ 1 → line plot per metric, groups overlaid (mean ± SEM)
        %
        % Inputs:
        %   sw      - struct from loadSlidingWindowData
        %   cfg     - PlotConfig
        %   filter  - GCMCResultsFilter
        %   metrics - cell array of metric column names (auto-detected if empty)
        %   opts    - optional name-value pairs:
        %       relative - logical, if true subtract pre-stimulus baseline (default: false)
        %       clim     - [nMetrics x 2] double, [cmin cmax] per metric row, or [] for auto
        %       average_by_subject - logical, if true average rows per subject (default: false)
        %
        % Outputs:
        %   hf         - figure handles: [nG x nM] (heatmap) or [nM x 1] (line)
        %   p_friedman - [nG x nW x nM] double, Friedman p per (group, window, metric)
        %               NaN where test is inapplicable

        arguments
            sw struct
            cfg PlotConfig = PlotConfig()
            filter GCMCResultsFilter = GCMCResultsFilter('shuffle', false)
            metrics cell = {}
            opts.relative logical = false
            opts.clim double = []
            opts.average_by_subject logical = false
        end

        nGroups = numel(sw.groups);
        nWindows = numel(sw.windows);

        % Auto-detect metrics from first non-empty table
        if isempty(metrics)
            for idx = 1:numel(sw.data)
                if ~isempty(sw.data{idx})
                    metrics = sw.data{idx}.Properties.VariableNames(5:end-3);
                    break;
                end
            end
        end
        nMetrics = numel(metrics);

        % Discover available reps across all filtered data
        reps = discoverReps(sw, filter);
        nReps = numel(reps);

        % Identify pre-stimulus windows (end time < 0)
        pre_stim_mask = sw.t_edges(:, 2) < 0;

        if nReps > 1
            [hf, p_friedman] = plotSWHeatmaps(sw, cfg, filter, metrics, reps, ...
                opts.relative, opts.clim, pre_stim_mask, opts.average_by_subject);
        else
            [hf, p_friedman] = plotSWLines(sw, cfg, filter, metrics, ...
                opts.relative, opts.clim, pre_stim_mask, opts.average_by_subject);
        end
    end

end

end

%% Local helper functions for plotMetricsByRepLine

function [data_by_group_rep, reps] = extractDataByGroupRep(group_data, filter, metric_name, average_by_subject)
    % Extract metric data organized by [group, rep]
    % If average_by_subject is true, return one value per subject (mean across rows)
    if nargin < 4; average_by_subject = false; end

    nGroups = numel(group_data);

    % Find available reps
    all_reps = [];
    for i = 1:nGroups
        filtered = filter.filterTable(group_data(i).data);
        if ismember('manifold_rep1', filtered.Properties.VariableNames)
            all_reps = unique([all_reps; filtered.manifold_rep1]);
        end
    end
    reps = sort(all_reps);
    nReps = numel(reps);

    % Extract data
    data_by_group_rep = cell(nGroups, nReps);
    for i_group = 1:nGroups
        for i_rep = 1:nReps
            rep_filter = filter;
            rep_filter.manifold_rep1 = reps(i_rep);
            rep_filter.manifold_rep2 = reps(i_rep);
            filtered = rep_filter.filterTable(group_data(i_group).data);
            if ~isempty(filtered) && ismember(metric_name, filtered.Properties.VariableNames)
                if average_by_subject
                    data_by_group_rep{i_group, i_rep} = averageMetricBySubject(filtered, metric_name);
                else
                    data_by_group_rep{i_group, i_rep} = filtered.(metric_name);
                end
            else
                data_by_group_rep{i_group, i_rep} = [];
            end
        end
    end
end

function [means, sems] = computeGroupRepStats(data_by_group_rep)
    % Compute mean and SEM for each cell in data_by_group_rep
    [nGroups, nReps] = size(data_by_group_rep);
    means = nan(nGroups, nReps);
    sems = nan(nGroups, nReps);

    for i = 1:nGroups
        for j = 1:nReps
            vals = data_by_group_rep{i, j};
            if ~isempty(vals)
                [means(i,j), sems(i,j)] = statsUtils.groupStats({vals});
            end
        end
    end
end

function friedman_results = computeFriedmanByGroup(group_data, filter, metric_name, reps)
    % Compute Friedman test for each group (subjects x reps)
    nGroups = numel(group_data);
    nReps = numel(reps);
    friedman_results = struct('group', {}, 'p', {}, 'chi2', {});

    for i_group = 1:nGroups
        friedman_results(i_group).group = group_data(i_group).group_name;
        friedman_results(i_group).p = NaN;
        friedman_results(i_group).chi2 = NaN;

        filtered = filter.filterTable(group_data(i_group).data);
        if ~ismember('subj_id', filtered.Properties.VariableNames)
            continue;
        end

        % Build subjects x reps matrix
        subj_ids = unique(filtered.subj_id);
        subj_rep_matrix = nan(numel(subj_ids), nReps);

        for i_subj = 1:numel(subj_ids)
            for i_rep = 1:nReps
                rep_filter = filter;
                rep_filter.manifold_rep1 = reps(i_rep);
                rep_filter.manifold_rep2 = reps(i_rep);
                rep_filter.subj_ids = subj_ids(i_subj);
                subj_data = rep_filter.filterTable(group_data(i_group).data);
                if ~isempty(subj_data) && ismember(metric_name, subj_data.Properties.VariableNames)
                    subj_rep_matrix(i_subj, i_rep) = mean(subj_data.(metric_name), 'omitnan');
                end
            end
        end

        % Run Friedman test
        result = statsUtils.friedman(subj_rep_matrix);
        friedman_results(i_group).p = result.p;
        friedman_results(i_group).chi2 = result.chi2;
    end
end

function between_results = computeBetweenGroupStats(data_by_group_rep, group_names, reps)
    % Mann-Whitney U between groups at each rep, with FDR correction
    [nGroups, nReps] = size(data_by_group_rep);
    between_results = struct('rep', {}, 'group1', {}, 'group2', {}, ...
                            'i_rep', {}, 'p_raw', {}, 'p_fdr', {});

    p_raw = [];
    idx = 0;

    for i_rep = 1:nReps
        groups_at_rep = data_by_group_rep(:, i_rep);
        mw_results = statsUtils.pairwiseMannWhitney(groups_at_rep, group_names, false);

        for k = 1:numel(mw_results)
            idx = idx + 1;
            between_results(idx).rep = reps(i_rep);
            between_results(idx).i_rep = i_rep;
            between_results(idx).group1 = mw_results(k).group1;
            between_results(idx).group2 = mw_results(k).group2;
            between_results(idx).p_raw = mw_results(k).p_raw;
            p_raw(idx) = mw_results(k).p_raw;
        end
    end

    % Apply FDR across all comparisons
    p_fdr = statsUtils.fdr(p_raw);
    for k = 1:numel(between_results)
        between_results(k).p_fdr = p_fdr(k);
    end
end

%% Local helper functions for sliding window methods

function [t_edges, valid] = parseSlidingWindowEdges(window_names)
    % Parse time window directory names into [t0, t1] edges
    %
    % Handles formats: '0.5-4.5s', '-4-0s', '-6--2s', '0_5-4_5s' (underscore for decimals)
    % Reuses regex pattern from sortByWindow in GCMC_Analysis.m
    %
    % Inputs:
    %   window_names - cell array of directory name strings
    %
    % Outputs:
    %   t_edges - [n x 2] double, [t0, t1] for each name
    %   valid   - [n x 1] logical, true if name was successfully parsed

    n = numel(window_names);
    t_edges = nan(n, 2);
    valid = false(n, 1);

    pattern = '([-]?\d+\.?\d*)-([-]?\d+\.?\d*)\w*';
    for i = 1:n
        name = strrep(window_names{i}, '_', '.'); % underscore→decimal convention
        tokens = regexp(name, pattern, 'tokens');
        if ~isempty(tokens)
            t0 = str2double(tokens{1}{1});
            t1 = str2double(tokens{1}{2});
            if ~isnan(t0) && ~isnan(t1)
                t_edges(i, :) = [t0, t1];
                valid(i) = true;
            end
        end
    end
end

function reps = discoverReps(sw, filter)
    % Find unique manifold_rep1 values across all filtered sliding window data
    reps = [];
    for i = 1:numel(sw.data)
        tbl = sw.data{i};
        if isempty(tbl); continue; end
        filtered = filter.filterTable(tbl);
        if ismember('manifold_rep1', filtered.Properties.VariableNames)
            reps = unique([reps; filtered.manifold_rep1]);
        end
    end
    reps = sort(reps);
end

function [hf, p_friedman] = plotSWHeatmaps(sw, cfg, filter, metrics, reps, do_relative, clim_arr, pre_stim_mask, average_by_subject)
    % Heatmap mode: one figure per (group, metric), dims [reps x windows]
    % Friedman test per (group, window): subjects × reps
    %
    % do_relative: if true, subtract pre-stimulus baseline per (subject, stim pair, rep)
    % clim_arr: [nMetrics x 2] double, or [] for auto per-figure
    % pre_stim_mask: logical [nWindows x 1], true for pre-stimulus windows
    % average_by_subject: if true, average by subject before computing heatmap cell means
    if nargin < 9; average_by_subject = false; end

    nGroups = numel(sw.groups);
    nWindows = numel(sw.windows);
    nMetrics = numel(metrics);
    nReps = numel(reps);

    hf = gobjects(nGroups, nMetrics);
    p_friedman = nan(nGroups, nWindows, nMetrics);

    for i_m = 1:nMetrics
        metric = metrics{i_m};

        % Get clim for this metric (if provided)
        if ~isempty(clim_arr) && size(clim_arr, 1) >= i_m
            metric_clim = clim_arr(i_m, :);
        else
            metric_clim = [];
        end

        for i_g = 1:nGroups
            if do_relative
                [heatmap_data, p_friedman(i_g, :, i_m)] = ...
                    computeHeatmapRelative(sw.data(i_g, :), filter, metric, reps, pre_stim_mask, average_by_subject);
            else
                [heatmap_data, p_friedman(i_g, :, i_m)] = ...
                    computeHeatmapAbsolute(sw.data(i_g, :), filter, metric, reps, average_by_subject);
            end

            % Plot
            hf(i_g, i_m) = figure;
            imagesc(sw.t_centers, reps, heatmap_data);
            if ~isempty(metric_clim)
                clim(metric_clim);
            end
            axis tight;
            xlabel('Time (s)');
            ylabel('Repetition #');
            suffix = '';
            if do_relative; suffix = ' (rel)'; end
            title(sprintf('%s - %s%s', sw.groups{i_g}, strrep(metric, 'Fun_', ''), suffix));
            colorbar;
            colormap(cfg.colormapName);
            set(gca, 'color', cfg.bgcol, 'XColor', cfg.axcol, 'YColor', cfg.axcol);
            set(gcf, 'color', cfg.bgcol);
            set(gcf, 'Position', [100, 100, 400, 300]);
        end
    end
end

function [hf, p_friedman] = plotSWLines(sw, cfg, filter, metrics, do_relative, clim_arr, pre_stim_mask, average_by_subject) %#ok<INUSD>
    % Line mode: one figure per metric, groups overlaid (mean ± SEM across subjects)
    % No Friedman test (no rep dimension)
    %
    % do_relative: if true, subtract pre-stimulus baseline per subject
    % clim_arr: [nMetrics x 2] double, or [] for auto
    % pre_stim_mask: logical [nWindows x 1], true for pre-stimulus windows
    % average_by_subject: accepted for API consistency (no behavior change;
    %   computeSubjectWindowMatrix already averages by subject)

    nGroups = numel(sw.groups);
    nWindows = numel(sw.windows);
    nMetrics = numel(metrics);

    hf = gobjects(nMetrics, 1);
    p_friedman = nan(nGroups, nWindows, nMetrics);

    for i_m = 1:nMetrics
        metric = metrics{i_m};

        % Get ylim for this metric (if provided)
        if ~isempty(clim_arr) && size(clim_arr, 1) >= i_m
            metric_ylim = clim_arr(i_m, :);
        else
            metric_ylim = [];
        end

        % Compute per-subject values for all groups
        mu_all = nan(nGroups, nWindows);
        sem_all = nan(nGroups, nWindows);

        for i_g = 1:nGroups
            subj_vals = computeSubjectWindowMatrix(sw.data(i_g, :), filter, metric);

            if do_relative && any(pre_stim_mask)
                baseline = mean(subj_vals(:, pre_stim_mask), 2, 'omitnan');
                subj_vals = subj_vals - baseline;
            end

            for i_w = 1:nWindows
                vals = subj_vals(:, i_w);
                [mu_all(i_g, i_w), sem_all(i_g, i_w)] = statsUtils.groupStats({vals});
            end
        end

        % Plot
        hf(i_m) = figure;
        hold on;
        line_handles = gobjects(nGroups, 1);

        for i_g = 1:nGroups
            mu = mu_all(i_g, :);
            sem = sem_all(i_g, :);
            color = cfg.c(i_g, :);
            valid = ~isnan(mu);

            if sum(valid) > 1
                lcfg.lineStyle = '-';
                lcfg.lineWidth = .5;
                plotLineNShade(sw.t_centers(valid), mu(valid)', sem(valid)', color, lcfg);
            end

            line_handles(i_g) = plot(sw.t_centers, mu, 'o-', ...
                'Color', color, 'MarkerFaceColor', color, ...
                'MarkerSize', 3, 'DisplayName', sw.groups{i_g});
        end

        hold off;
        xlabel('Time (s)');
        suffix = '';
        if do_relative; suffix = ' (rel)'; end
        ylabel([strrep(metric, 'Fun_', ''), suffix]);
        if ~isempty(metric_ylim)
            ylim(metric_ylim);
        end
        legend(line_handles, 'Location', 'best', 'Box', 'off');
        box off;
        set(gca, 'color', cfg.bgcol, 'XColor', cfg.axcol, 'YColor', cfg.axcol);
        set(gcf, 'color', cfg.bgcol);
        set(gcf, 'Position', [100, 100, 400, 350]);
    end
end

function p = friedmanAtWindow(filtered, metric, reps)
    % Friedman test at a single window: subjects × reps
    p = NaN;
    nReps = numel(reps);
    if ~ismember('subj_id', filtered.Properties.VariableNames) || nReps < 2
        return;
    end

    subj_ids = unique(filtered.subj_id);
    nSubj = numel(subj_ids);
    if nSubj < 2; return; end

    subj_rep_mat = nan(nSubj, nReps);
    for i_s = 1:nSubj
        for i_r = 1:nReps
            mask = filtered.subj_id == subj_ids(i_s) & ...
                   filtered.manifold_rep1 == reps(i_r) & ...
                   filtered.manifold_rep2 == reps(i_r);
            vals = filtered{mask, metric};
            if ~isempty(vals)
                subj_rep_mat(i_s, i_r) = mean(vals, 'omitnan');
            end
        end
    end

    result = statsUtils.friedman(subj_rep_mat);
    p = result.p;
end

function [heatmap_data, p_friedman] = computeHeatmapAbsolute(data_row, filter, metric, reps, average_by_subject)
    % Compute heatmap data (absolute values) for a single group
    % data_row: {1 x nWindows} cell array of tables
    % average_by_subject: if true, average by subject before taking cell mean
    % Returns: heatmap_data [nReps x nWindows], p_friedman [1 x nWindows]
    if nargin < 5; average_by_subject = false; end

    nWindows = numel(data_row);
    nReps = numel(reps);
    heatmap_data = nan(nReps, nWindows);
    p_friedman = nan(1, nWindows);

    for i_w = 1:nWindows
        tbl = data_row{i_w};
        if isempty(tbl); continue; end
        filtered = filter.filterTable(tbl);
        if isempty(filtered) || ~ismember(metric, filtered.Properties.VariableNames)
            continue;
        end

        for i_r = 1:nReps
            mask = filtered.manifold_rep1 == reps(i_r) & ...
                   filtered.manifold_rep2 == reps(i_r);
            rep_rows = filtered(mask, :);
            if ~isempty(rep_rows)
                if average_by_subject
                    vals = averageMetricBySubject(rep_rows, metric);
                else
                    vals = rep_rows.(metric);
                end
                heatmap_data(i_r, i_w) = mean(vals, 'omitnan');
            end
        end

        p_friedman(i_w) = friedmanAtWindow(filtered, metric, reps);
    end
end

function [heatmap_data, p_friedman] = computeHeatmapRelative(data_row, filter, metric, reps, pre_stim_mask, average_by_subject)
    % Compute heatmap data relative to pre-stimulus baseline
    % Baseline subtraction happens per (subject, stim-pair, rep) BEFORE averaging
    %
    % data_row: {1 x nWindows} cell array of tables
    % pre_stim_mask: logical [nWindows x 1], true for pre-stimulus windows
    % average_by_subject: if true, average by subject after baseline subtraction
    % Returns: heatmap_data [nReps x nWindows], p_friedman [1 x nWindows]
    if nargin < 6; average_by_subject = false; end

    nWindows = numel(data_row);
    nReps = numel(reps);
    heatmap_data = nan(nReps, nWindows);
    p_friedman = nan(1, nWindows);

    if ~any(pre_stim_mask)
        % No pre-stimulus windows → fall back to absolute
        [heatmap_data, p_friedman] = computeHeatmapAbsolute(data_row, filter, metric, reps, average_by_subject);
        return;
    end

    % Build a map: (subj_id, stim_pair_key, rep) → baseline value
    % stim_pair_key = sprintf('%d_%d', manifold_idx_1, manifold_idx_2)
    baseline_map = containers.Map('KeyType', 'char', 'ValueType', 'double');

    % Collect baseline values from pre-stimulus windows
    pre_win_idx = find(pre_stim_mask);
    for i_w = pre_win_idx(:)'
        tbl = data_row{i_w};
        if isempty(tbl); continue; end
        filtered = filter.filterTable(tbl);
        if isempty(filtered) || ~ismember(metric, filtered.Properties.VariableNames)
            continue;
        end
        if ~ismember('subj_id', filtered.Properties.VariableNames); continue; end

        for i_row = 1:height(filtered)
            row = filtered(i_row, :);
            key = sprintf('%d_%d_%d_%d', row.subj_id, row.manifold_idx_1, ...
                row.manifold_idx_2, row.manifold_rep1);
            val = row.(metric);
            if isKey(baseline_map, key)
                baseline_map(key) = [baseline_map(key), val];
            else
                baseline_map(key) = val;
            end
        end
    end

    % Average baselines
    keys = baseline_map.keys;
    for i = 1:numel(keys)
        baseline_map(keys{i}) = mean(baseline_map(keys{i}), 'omitnan');
    end

    % Compute relative values for all windows
    for i_w = 1:nWindows
        tbl = data_row{i_w};
        if isempty(tbl); continue; end
        filtered = filter.filterTable(tbl);
        if isempty(filtered) || ~ismember(metric, filtered.Properties.VariableNames)
            continue;
        end
        if ~ismember('subj_id', filtered.Properties.VariableNames); continue; end

        for i_r = 1:nReps
            mask = filtered.manifold_rep1 == reps(i_r) & ...
                   filtered.manifold_rep2 == reps(i_r);
            rows = filtered(mask, :);
            if isempty(rows); continue; end

            rel_vals = nan(height(rows), 1);
            for j = 1:height(rows)
                row = rows(j, :);
                key = sprintf('%d_%d_%d_%d', row.subj_id, row.manifold_idx_1, ...
                    row.manifold_idx_2, row.manifold_rep1);
                if isKey(baseline_map, key)
                    rel_vals(j) = row.(metric) - baseline_map(key);
                else
                    rel_vals(j) = row.(metric);  % no baseline available
                end
            end

            if average_by_subject && ismember('subj_id', rows.Properties.VariableNames)
                % Average relative values per subject, then take mean
                subj_ids = unique(rows.subj_id);
                subj_means = nan(numel(subj_ids), 1);
                for js = 1:numel(subj_ids)
                    subj_means(js) = mean(rel_vals(rows.subj_id == subj_ids(js)), 'omitnan');
                end
                heatmap_data(i_r, i_w) = mean(subj_means, 'omitnan');
            else
                heatmap_data(i_r, i_w) = mean(rel_vals, 'omitnan');
            end
        end

        % Friedman on relative values requires re-filtering with baseline subtraction
        % For simplicity, use Friedman on absolute (test is about rank differences anyway)
        p_friedman(i_w) = friedmanAtWindow(filtered, metric, reps);
    end
end

function subj_window_mat = computeSubjectWindowMatrix(data_row, filter, metric)
    % Compute [nSubjects x nWindows] matrix of per-subject mean values
    % data_row: {1 x nWindows} cell array of tables

    nWindows = numel(data_row);

    % First pass: find all subject IDs
    all_subj = [];
    for i_w = 1:nWindows
        tbl = data_row{i_w};
        if isempty(tbl); continue; end
        filtered = filter.filterTable(tbl);
        if ~isempty(filtered) && ismember('subj_id', filtered.Properties.VariableNames)
            all_subj = unique([all_subj; filtered.subj_id]);
        end
    end
    nSubj = numel(all_subj);
    subj_window_mat = nan(nSubj, nWindows);

    % Second pass: fill matrix
    for i_w = 1:nWindows
        tbl = data_row{i_w};
        if isempty(tbl); continue; end
        filtered = filter.filterTable(tbl);
        if isempty(filtered) || ~ismember(metric, filtered.Properties.VariableNames)
            continue;
        end

        if ismember('subj_id', filtered.Properties.VariableNames)
            for i_s = 1:nSubj
                sv = filtered{filtered.subj_id == all_subj(i_s), metric};
                subj_window_mat(i_s, i_w) = mean(sv, 'omitnan');
            end
        else
            % No subj_id column → treat all rows as one "subject"
            subj_window_mat(1, i_w) = mean(filtered.(metric), 'omitnan');
        end
    end
end

function vals = averageMetricBySubject(filtered, metric_name)
    % Average metric values per unique subject, returning one value per subject
    %
    % Inputs:
    %   filtered    - table with rows to average
    %   metric_name - name of the metric column
    %
    % Output:
    %   vals - [nSubjects x 1] vector of per-subject means
    %          Falls back to raw values if no subj_id column exists

    if ~ismember('subj_id', filtered.Properties.VariableNames)
        vals = filtered.(metric_name);
        return;
    end

    subj_ids = unique(filtered.subj_id);
    nSubj = numel(subj_ids);
    vals = nan(nSubj, 1);
    for i = 1:nSubj
        vals(i) = mean(filtered{filtered.subj_id == subj_ids(i), metric_name}, 'omitnan');
    end
end