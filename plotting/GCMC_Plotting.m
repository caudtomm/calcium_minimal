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

    function 

    function plotBoxplotsByGroup(group_data, cfg, shuffle)
        arguments
            group_data struct
            cfg PlotConfig = PlotConfig()
            shuffle logical = false;
        end

        nGroups = numel(group_data);
        metrics = group_data(1).data.Properties.VariableNames(5:end-3); % exclude grouping variables and stimulus names
        nMetrics = numel(metrics);

        for i_metric = 1:7%nMetrics
            figure;
            hold on;

            disp(['Metric: ', metrics{i_metric}])

            % Collect data for each group
            group_labels = {};
            box_data = [];
            datacells = cell(1,nGroups);
            for i_group = 1:nGroups
                this_data = group_data(i_group).data{shuffle==group_data(i_group).data.shuffle, metrics{i_metric}};
                box_data = [box_data; this_data];
                datacells{i_group} = this_data(:);
                group_labels = [group_labels; repelem(string(group_data(i_group).group_name), size(this_data, 1), 1)];
            end

            % Mann-Whitney U-test (non-parametric test)
            for i_group = 1:nGroups
                this_data = group_data(i_group).data{shuffle==group_data(i_group).data.shuffle, metrics{i_metric}};
                for j_group = i_group+1:nGroups
                    other_data = group_data(j_group).data{shuffle==group_data(j_group).data.shuffle, metrics{i_metric}};
                    p = ranksum(this_data, other_data); % Mann-Whitney U-test
                    disp(['Mann-Whitney U-test between ', group_data(i_group).group_name, ...
                            ' and ', group_data(j_group).group_name, ...
                            ' for metric ', metrics{i_metric}, ': p = ', num2str(p)]);
                end
            end

            % Create boxplot
            %boxplot(box_data, group_labels, 'Notch', 'on', 'Labels', unique(group_labels, 'stable'));
            RF_mkBoxPlot3(datacells,[],[],.5,1,2,15,[]);
            xticks([1:nGroups]); xticklabels({group_data(:).group_name})
            % Superimpose scatter plot for each group
            for i_group = 1:nGroups
                this_data = group_data(i_group).data{shuffle==group_data(i_group).data.shuffle, metrics{i_metric}};
                this_subj_ids = group_data(i_group).data{shuffle==group_data(i_group).data.shuffle,'subj_id'};
                scatter(repelem(i_group, numel(this_data)), ...
                        this_data, 40, 'filled', 'CData', cfg.c(this_subj_ids,:), 'MarkerFaceAlpha', 0.7, ...
                        'jitter', 'on', 'jitterAmount', 0.15);
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

    function plotBoxplotsForEachMetric(avg_results, cfg)
        arguments
            avg_results table % from a single subj
            cfg PlotConfig = PlotConfig()
        end

        metrics = avg_results.Properties.VariableNames(5:end-2); % exclude grouping variables and stimulus names

        % separate data and shuffle conditions
        y_data = avg_results{avg_results.shuffle==false, metrics};
        y_shuffle = avg_results{avg_results.shuffle==true, metrics};

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

    function plotMetricStability(results, cfg)
        arguments
            results table
            cfg PlotConfig = PlotConfig()
        end


        results = results(results.shuffle==false, :); % only data, no shuffle

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

end

end