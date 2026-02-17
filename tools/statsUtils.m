classdef statsUtils
    % statsUtils - Statistical utility functions
    %
    % Static methods for common statistical operations:
    %   - p-value formatting
    %   - multiple comparison corrections
    %   - non-parametric tests with structured output
    %
    % Example usage:
    %   stars = statsUtils.pvalToStars(0.003);  % '**'
    %   p_adj = statsUtils.fdr([0.01, 0.04, 0.06]);  % FDR correction
    %   result = statsUtils.friedman(data_matrix);  % Friedman test

    methods (Static)

        function stars = pvalToStars(p, thresholds)
            % Convert p-value to significance stars
            %
            % Inputs:
            %   p - p-value (scalar or array)
            %   thresholds - [optional] struct with fields p001, p01, p05
            %                default: 0.001, 0.01, 0.05
            %
            % Output:
            %   stars - string: '***', '**', '*', or 'n.s.'

            arguments
                p double
                thresholds.p001 double = 0.001
                thresholds.p01 double = 0.01
                thresholds.p05 double = 0.05
            end

            if isscalar(p)
                if isnan(p)
                    stars = '';
                elseif p < thresholds.p001
                    stars = '***';
                elseif p < thresholds.p01
                    stars = '**';
                elseif p < thresholds.p05
                    stars = '*';
                else
                    stars = 'n.s.';
                end
            else
                % Handle array input
                stars = cell(size(p));
                for i = 1:numel(p)
                    stars{i} = statsUtils.pvalToStars(p(i), ...
                        'p001', thresholds.p001, 'p01', thresholds.p01, 'p05', thresholds.p05);
                end
            end
        end

        function p_adj = fdr(p_vals)
            % Benjamini-Hochberg FDR correction
            %
            % Input:
            %   p_vals - vector of raw p-values
            %
            % Output:
            %   p_adj - FDR-adjusted p-values (same size as input)

            p_adj = nan(size(p_vals));
            valid = ~isnan(p_vals);

            if sum(valid) == 0
                return;
            end

            p_valid = p_vals(valid);
            n = numel(p_valid);
            [p_sorted, sort_idx] = sort(p_valid);

            % BH procedure: p_adj(i) = p(i) * n / rank(i)
            ranks = (1:n)';
            adjusted = p_sorted(:) .* n ./ ranks;

            % Enforce monotonicity (cumulative minimum from end)
            for i = n-1:-1:1
                adjusted(i) = min(adjusted(i), adjusted(i+1));
            end
            adjusted = min(adjusted, 1);  % Cap at 1

            % Restore original order
            p_adj_unsorted = nan(size(p_valid));
            p_adj_unsorted(sort_idx) = adjusted;
            p_adj(valid) = p_adj_unsorted;
        end

        function result = friedman(data_matrix, display_opt)
            % Friedman test with structured output
            %
            % Input:
            %   data_matrix - [subjects x conditions] matrix
            %   display_opt - 'on' or 'off' (default: 'off')
            %
            % Output:
            %   result - struct with fields:
            %       .p      - p-value
            %       .chi2   - chi-square statistic
            %       .df     - degrees of freedom
            %       .n      - number of subjects
            %       .k      - number of conditions
            %       .valid  - true if test was performed

            arguments
                data_matrix double
                display_opt char {mustBeMember(display_opt, {'on', 'off'})} = 'off'
            end

            result = struct('p', NaN, 'chi2', NaN, 'df', NaN, ...
                           'n', 0, 'k', 0, 'valid', false);

            % Remove rows with NaN
            valid_rows = all(~isnan(data_matrix), 2);
            data_clean = data_matrix(valid_rows, :);

            [n, k] = size(data_clean);
            result.n = n;
            result.k = k;

            if n < 2 || k < 2
                return;
            end

            try
                [p, tbl, stats] = friedman(data_clean, 1, display_opt);
                result.p = p;
                % result.chi2 = stats.chisq;
                result.df = k - 1;
                result.valid = true;
            catch
                % Test failed
            end
        end

        function results = pairwiseMannWhitney(groups, group_names, apply_fdr)
            % Pairwise Mann-Whitney U tests between groups
            %
            % Inputs:
            %   groups - cell array of data vectors, one per group
            %   group_names - cell array of group name strings
            %   apply_fdr - logical, whether to apply FDR correction (default: true)
            %
            % Output:
            %   results - struct array with fields:
            %       .group1, .group2 - group names
            %       .i, .j - group indices
            %       .p_raw - raw p-value
            %       .p_adj - FDR-adjusted p-value (if apply_fdr)
            %       .significant - true if p_adj < 0.05

            arguments
                groups cell
                group_names cell = {}
                apply_fdr logical = true
            end

            nGroups = numel(groups);
            if isempty(group_names)
                group_names = arrayfun(@(x) sprintf('Group%d', x), 1:nGroups, 'UniformOutput', false);
            end

            % Count number of comparisons
            nComparisons = nGroups * (nGroups - 1) / 2;
            results = struct('group1', {}, 'group2', {}, 'i', {}, 'j', {}, ...
                           'p_raw', {}, 'p_adj', {}, 'significant', {});

            p_raw = nan(nComparisons, 1);
            comp_idx = 0;

            for i = 1:nGroups
                for j = (i+1):nGroups
                    comp_idx = comp_idx + 1;

                    data1 = groups{i};
                    data2 = groups{j};

                    if numel(data1) >= 2 && numel(data2) >= 2
                        p_raw(comp_idx) = ranksum(data1(:), data2(:));
                    else
                        p_raw(comp_idx) = NaN;
                    end

                    results(comp_idx).group1 = group_names{i};
                    results(comp_idx).group2 = group_names{j};
                    results(comp_idx).i = i;
                    results(comp_idx).j = j;
                    results(comp_idx).p_raw = p_raw(comp_idx);
                end
            end

            % Apply FDR correction
            if apply_fdr
                p_adj = statsUtils.fdr(p_raw);
            else
                p_adj = p_raw;
            end

            for k = 1:nComparisons
                results(k).p_adj = p_adj(k);
                results(k).significant = p_adj(k) < 0.05;
            end
        end

        function [means, sems, ns] = groupStats(data_cell)
            % Compute mean and SEM for each group
            %
            % Input:
            %   data_cell - cell array of data vectors
            %
            % Output:
            %   means - vector of means
            %   sems - vector of standard errors
            %   ns - vector of sample sizes

            nGroups = numel(data_cell);
            means = nan(nGroups, 1);
            sems = nan(nGroups, 1);
            ns = nan(nGroups, 1);

            for i = 1:nGroups
                vals = data_cell{i};
                vals = vals(~isnan(vals));
                ns(i) = numel(vals);
                if ns(i) > 0
                    means(i) = mean(vals);
                    sems(i) = std(vals) / sqrt(ns(i));
                end
            end
        end

    end
end
