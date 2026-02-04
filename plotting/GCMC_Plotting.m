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

    function hf = plotBoxplotsByGroup(group_data, cfg, filter)
        arguments
            group_data struct
            cfg PlotConfig = PlotConfig()
            filter GCMCResultsFilter = GCMCResultsFilter('shuffle', false)
        end

        nGroups = numel(group_data);
        metrics = group_data(1).data.Properties.VariableNames(5:end-3); % exclude grouping variables and stimulus names
        nMetrics = 7; %numel(metrics);

        hf = gobjects(nMetrics,1);

        for i_metric = 1:nMetrics
            hf(i_metric) = figure;
            hold on;

            disp(['Metric: ', metrics{i_metric}])

            % Collect data for each group
            group_labels = {};
            box_data = [];
            datacells = cell(1,nGroups);
            for i_group = 1:nGroups
                filtered_table = filter.filterTable(group_data(i_group).data);
                this_data = filtered_table{:, metrics{i_metric}};
                box_data = [box_data; this_data];
                datacells{i_group} = this_data(:);
                group_labels = [group_labels; repelem(string(group_data(i_group).group_name), size(this_data, 1), 1)];
            end

            % Mann-Whitney U-test (non-parametric test)
            for i_group = 1:nGroups
                filtered_table = filter.filterTable(group_data(i_group).data);
                this_data = filtered_table{:, metrics{i_metric}};
                for j_group = i_group+1:nGroups
                    filtered_other = filter.filterTable(group_data(j_group).data);
                    other_data = filtered_other{:, metrics{i_metric}};
                    p = ranksum(this_data, other_data); % Mann-Whitney U-test
                    disp(['Mann-Whitney U-test between ', group_data(i_group).group_name, ...
                            ' and ', group_data(j_group).group_name, ...
                            ' for metric ', metrics{i_metric}, ': p = ', num2str(p)]);
                end
            end

            % Create boxplot
            %boxplot(box_data, group_labels, 'Notch', 'on', 'Labels', unique(group_labels, 'stable'));
            RF_mkBoxPlot3(datacells,[],[],.1,.5,.5,2,[]);
            xticks([1:nGroups]); xticklabels({group_data(:).group_name})
            % Superimpose scatter plot for each group
            for i_group = 1:nGroups
                filtered_table = filter.filterTable(group_data(i_group).data);
                this_data = filtered_table{:, metrics{i_metric}};
                this_subj_ids = filtered_table{:, 'subj_id'};
                %scatter(repelem(i_group, numel(this_data)), ...
                %        this_data, 40, 'filled', 'CData', cfg.c(this_subj_ids,:), 'MarkerFaceAlpha', 0.7, ...
                %        'jitter', 'on', 'jitterAmount', 0.15);
            end
            title(['Metric: ', metrics{i_metric}]);
            ylabel(metrics{i_metric});

            % Customize plot appearance
            box off;
            axis tight
            xlim([.5 max(xticks)+.5])
            set(gca, 'color', cfg.bgcol, 'XColor', cfg.axcol, 'YColor', cfg.axcol, 'ZColor', cfg.axcol);
            set(gcf, 'color', cfg.bgcol);
            set(gcf, 'Position', [100, 100, 200, 500]);
            hold off;
        end

    end

    function plotBoxplotsForEachMetric(avg_results, cfg, filter)
        arguments
            avg_results table % from a single subj
            cfg PlotConfig = PlotConfig()
            filter GCMCResultsFilter = GCMCResultsFilter()
        end

        metrics = avg_results.Properties.VariableNames(5:end-2); % exclude grouping variables and stimulus names

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

            thisdata = y_data(:,i);
            thisshuffle = y_shuffle(:,i);

            % stats

            % determine which tail to use based on difference in means.
            % this makes sense only because the shuffle is the internal control!
            % it's just a quicker way to define what's expected to be higher or lower in the data.
            mean_data = mean(thisdata);
            mean_shuffle = mean(thisshuffle);
            if mean_data > mean_shuffle
                tail = 'right';
            else
                tail = 'left';
            end

            try
                [~, p] = ttest(thisdata, thisshuffle, 'Tail', tail);
            catch
                p = NaN;
            end
            disp(['Metric: ', metrics{i}, ', p-value (one-tailed paired t-test): ', num2str(p)]);

            boxplot([thisdata;thisshuffle], ...
                [repelem("Data", size(y_data,1),1); repelem("Shuffle", size(y_shuffle,1),1)]);
            hold on;
            scatter(repelem(1, size(y_data, 1)), y_data(:, i), 'r', 'filled', 'jitter', 'on', 'jitterAmount', 0.15);
            scatter(repelem(2, size(y_shuffle, 1)), y_shuffle(:, i), 'b', 'filled', 'jitter', 'on', 'jitterAmount', 0.15);
            hold off;
            title([metrics{i}, ', p=', sprintf(num2str(p), '%.5f')]);
            ylabel(metrics{i});

            box off
            set(gcf,'Position',[100 100 200 400])
            set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
            set(gcf, 'color', cfg.bgcol); 
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

    function [hf, stats] = plotMetricsByRepLine(group_data, cfg, filter, metric_name)
        % Plot metric values across repetitions as line plot with error ribbons
        %
        % Uses statsUtils for statistics and plotLineNShade for visualization.
        %
        % Inputs:
        %   group_data - struct array with .group_name and .data (table)
        %   cfg - PlotConfig object
        %   filter - GCMCResultsFilter
        %   metric_name - name of metric column (e.g., 'Fun_capacity')
        %
        % Outputs:
        %   hf - figure handle
        %   stats - struct with .friedman and .between_group results

        arguments
            group_data struct
            cfg PlotConfig = PlotConfig()
            filter GCMCResultsFilter = GCMCResultsFilter('shuffle', false)
            metric_name char = 'Fun_capacity'
        end

        nGroups = numel(group_data);
        group_names = {group_data.group_name};

        % === Extract data by group and rep ===
        [data_by_group_rep, reps] = extractDataByGroupRep(group_data, filter, metric_name);
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
                lcfg.lineWidth = 2;
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
                'MarkerSize', 6, 'DisplayName', legend_str);
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

end

end

%% Local helper functions for plotMetricsByRepLine

function [data_by_group_rep, reps] = extractDataByGroupRep(group_data, filter, metric_name)
    % Extract metric data organized by [group, rep]
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
                data_by_group_rep{i_group, i_rep} = filtered.(metric_name);
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