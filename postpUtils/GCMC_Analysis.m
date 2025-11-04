classdef GCMC_Analysis
    properties
        viewer ExperimentViewer
    end
    methods
        function obj = GCMC_Analysis(viewer)
            arguments
                viewer ExperimentViewer
            end
            obj.viewer = viewer;
        end

        function manifolds = outputDataFiles(obj, outdir) % assumption: manifolds are 1 per stimulus per subject
            arguments
                obj GCMC_Analysis
                outdir char = 'manifold_data';
            end
            % prepare data for GCMC analysis
            % single fish
            % 
            v = obj.viewer;

            %% getting a list of the stimuli
            [~,labs] = v.dataFilter.filterData(v);
            stims = unique(labs{1});
            nStims = numel(stims);
            nSubjects = numel(labs);

            %% data extraction
            disp('Extracting manifolds for GCMC analysis...');
            manifolds = cell(nSubjects,nStims);
            for i_stim = 1:nStims
                % only get this stimulus's responses
                v.dataFilter.stims_allowed = stims(i_stim);
                thispoints = v.dataFilter.filterData(v);

                % concatenate all repetitions (one stimulus = one manifold)
                thispoints = cellfun(@ActivityTraces.format, thispoints, 'UniformOutput',false);

                % eliminate any NaN
                idx_tokeep = cellfun(@(x) sum(isnan(x),2)==0, thispoints, 'UniformOutput',false);
                for i_sub = 1:nSubjects
                    thispoints{i_sub} = thispoints{i_sub}(idx_tokeep{i_sub},:);
                end

                % [t,N] -> [N,t]
                thispoints = cellfun(@transpose, thispoints, 'UniformOutput',false);

                % store to manifolds
                manifolds(:,i_stim) = thispoints;
            end

            %% save each subject's data to a python-compatible mat file
            disp(['Saving manifold data to ', outdir, '...']);
            if ~isfolder(outdir); mkdir(outdir); end
            allSubjIDs = v.dataFilter.getSubjectIDs(v.subjectTab);
            for i_sub = 1:nSubjects
                data = manifolds(i_sub,:);
                subjID = allSubjIDs{i_sub};
                save(fullfiletol(outdir,['manifolds_subj',num2str(i_sub),'.mat']), "data", "stims", "subjID");
            end

            disp('Done.');

        end

        function [results, avg_results] = extractResults(obj,indir)
            arguments
                obj GCMC_Analysis
                indir char = pwd
            end

            v = obj.viewer;

            %% prepare output file
            outfname = fullfiletol(indir, 'GCMC_results_table.mat');
            if isfile(outfname)
                disp(['Existing GCMC results table found: ', outfname,'. Overwriting...']);
                delete(outfname);
            end

            %% read and parse results

            % getting stimulus names ( # TODO : this should be retrieved from the output files,
            % which in turn should inherit it from their input file
            v.dataFilter.subjectIDs = {'TC_240104_TC0028_240101beh1b3_sxpDp_odorexp004_RPB3144501500AG'};
            v.dataFilter.stims_allowed = 'all stimuli';
            [~,labs] = v.dataFilter.filterData(v);
            stims = unique(labs{1});

            files = dir(fullfiletol(indir,'*.mat'));

            % extract all results from individual files (pairwise manifold comparisons)
            allresults = cell(numel(files),1);
            for i = 1:numel(files)
                disp(['loading ',files(i).name])

                p = parseFileName(files(i).name);
                data = load(fullfiletol(indir,files(i).name));
                
                data.manifold_idx_1 = p.manifold_idx_1;
                data.manifold_idx_2 = p.manifold_idx_2;
                allresults{i} = data;
            end

            %% combine into a table
            results = table;
            for i = 1:numel(allresults)
                t = parseMetrics2table(allresults{i});
                
                % add the stimulus names
                t.manifold_name_1 = stims(t.manifold_idx_1 +1);
                t.manifold_name_2 = stims(t.manifold_idx_2 +1);

                results = [results; t];
            end

            %% average over repetitions with the same parameters
            avg_results = genResultAvgTable(results, [], stims);

            %% store results
            save(outfname, 'results', 'avg_results');
            disp(['Saved GCMC results table to ', outfname]);
        end

        function [all_results] = extractResultsFromMultipleSubjects(obj, indir)
            arguments
                obj GCMC_Analysis
                indir char = pwd
            end

            files = dir(fullfiletol(indir,'manifolds_subj*.mat'));
            nSubjects = numel(files);

            all_results = struct();
            counter = 0;
            for i_sub = 1:nSubjects
                subj_name = ['manifolds_subj',num2str(i_sub)];
                subdir = fullfiletol(indir,subj_name);
                if exist(subdir,'dir')==0
                    disp(['Skipping subject ', num2str(i_sub), ': folder not found - ', subdir]);
                    continue;
                end
                counter = counter+1;
                disp(['Extracting results for subject ', num2str(i_sub), ' from ', subdir]);

                [results, avg_results] = obj.extractResults(subdir);
                all_results(counter).name = subj_name;
                all_results(counter).num_id = i_sub;
                all_results(counter).results = results;
                all_results(counter).avg_results = avg_results;

            end
        end

        function plotInterGroupComparison(obj, all_results, cfg)
            arguments
                obj GCMC_Analysis
                all_results struct
                cfg PlotConfig = PlotConfig()
            end
            
            metrics = all_results(1).avg_results.Properties.VariableNames(5:end-2); % exclude grouping variables and stimulus names
            subjectTab = obj.viewer.subjectTab;
            nResults = numel(all_results);
            
            % select subjectTab rows (subjects) according to the results'
            % ordinal ID
            res_ids = [];
            for i = 1:nResults
                this_id = all_results(i).num_id;
                res_ids = [res_ids; this_id];
            end
            subjectTab = subjectTab(res_ids,:);

            % separate data by subject group
            groups = unique(subjectTab.group);
            nGroups = numel(groups);
            group_data = struct();
            for i = 1:nGroups
                idx = ismember(subjectTab.group,groups{i});
                this_res = all_results(idx);
                group_data(i).group_name = groups{i};
                group_data(i).data = table();

                for i_subj = 1:numel(this_res)
                    this_data = this_res(i_subj).avg_results;
                    group_data(i).data = [group_data(i).data; this_data]
                end

            end


        end

    end

    methods (Static)
        %% plotting

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

        function plotBoxplotsForEachMetric(avg_results, cfg)
            arguments
                avg_results table
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

%% helper functions


function avg_results = genResultAvgTable(results, idx, stims)
    if nargin < 2 || isempty(idx)
        avg_results = varfun(@(x) mean(x,"omitmissing"), results, ...
            'InputVariables', results.Properties.VariableNames(6:end-2), ...
            'GroupingVariables', {'manifold_idx_1','manifold_idx_2','shuffle'});
    else
        avg_results = varfun(@(x) mean(x(idx,:),"omitmissing"), results, ...
            'InputVariables', results.Properties.VariableNames(6:end-2), ...
            'GroupingVariables', {'manifold_idx_1','manifold_idx_2','shuffle'});
    end

    % add the stimulus names
    if ~exist("stims","var") || isempty(stims); return; end
    avg_results.manifold_name_1 = stims(avg_results.manifold_idx_1 +1);
    avg_results.manifold_name_2 = stims(avg_results.manifold_idx_2 +1);
end

function T = parseMetrics2table(s)
    % s : results structure

    M = s.metrics;

    % front-append index names tp column names
    cols = string(M.columns);     % 1×K string array (column labels)
    idx_names = {'manifold_idx_1','manifold_idx_2','repetition','shuffle','seed'};
    cols = [idx_names, cols];

    % front-append index values to data matrix
    idx  = M.index;               % Nx1 (char cell) - format {'(rep,shuffle,seed)'}
    [rep,shuf,seed] = cellfun(@(x) parseIndexVal(x), idx, 'UniformOutput', false);
    rep = cell2mat(rep); shuf = cell2mat(shuf); seed = cell2mat(seed);
    N = numel(rep);
    mid1 = repelem(s.manifold_idx_1,N,1); mid2 = repelem(s.manifold_idx_2,N,1);
    X    = M.values;              % NxK double
    X    = [mid1,mid2,rep,shuf,seed, X]; % Nx(K+7) double
    
    T = array2table(X, 'VariableNames', cellstr(cols));
    
end

function [firstInt, boolVal, lastInt] = parseIndexVal(s)
    % s : char vector or string
    
    % Extract integers
    nums = regexp(s, '\d+', 'match');
    firstInt = str2double(nums{1});
    lastInt  = str2double(nums{end});
    
    % Extract boolean
    boolStr = regexp(s, '(True|False)', 'match', 'ignorecase');
    boolVal = strcmpi(boolStr{1}, 'True');
end

function out = parseFileName(fname)
    % expected format: p_<manifold_idx_1>_<manifold_idx_2>.mat
    tokens = regexp(fname, '^p_(\d+)_(\d+)\.mat$', 'tokens');
    if isempty(tokens)
        error('Filename does not match expected format: %s', fname);
    end

    % extract indices
    tokens = tokens{1};
    out.manifold_idx_1 = str2double(tokens{1});
    out.manifold_idx_2 = str2double(tokens{2});
end