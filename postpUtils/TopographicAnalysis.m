classdef TopographicAnalysis
    properties
        at          % ActivityTraces object
        centroids   % (N, 2) centroid positions [x, y] for good neurons
        roi_ids     % goodNeuron_IDs used
    end

    methods
        function obj = TopographicAnalysis(activityTraces)
            arguments
                activityTraces ActivityTraces
            end
            obj.at = activityTraces;
            obj.roi_ids = activityTraces.goodNeuron_IDs;
            obj = obj.getCentroids();
        end

        function obj = getCentroids(obj)
            map = obj.at.ROImap;
            ids = obj.roi_ids;
            N = numel(ids);
            obj.centroids = zeros(N, 2);
            for i = 1:N
                props = regionprops(double(map == ids(i)), 'Centroid');
                obj.centroids(i, :) = props.Centroid; % [x, y]
            end
        end

        function ax = plotTopoMap(obj, unit_metric, ax)
            arguments
                obj
                unit_metric (:,1) double
                ax = []
            end
            assert(numel(unit_metric) == numel(obj.roi_ids), ...
                'unit_metric length (%d) must match number of good neurons (%d)', ...
                numel(unit_metric), numel(obj.roi_ids));

            map = obj.at.ROImap;
            img = nan(size(map));
            for i = 1:numel(obj.roi_ids)
                img(map == obj.roi_ids(i)) = unit_metric(i);
            end

            if isempty(ax); figure; ax = gca; end
            imagesc(ax, img, 'AlphaData', ~isnan(img));
            set(ax, 'Color', [0 0 0]);
            axis(ax, 'image'); axis(ax, 'off');
            colormap(ax, 'parula');
            colorbar(ax);
        end

        function [hf, tuningCurves, quant] = plotFavouriteStimulus(obj, events, stim_types, cfg, nShuffles, binSize, doPlot)
        % PLOTFAVOURITESTIMULUS  ROI maps coloured by favourite stimulus + distance quantification.
        %
        %   [hf, tuningCurves, quant] = obj.plotFavouriteStimulus(events, stim_types)
        %   [hf, tuningCurves, quant] = obj.plotFavouriteStimulus(events, stim_types, cfg)
        %   [hf, tuningCurves, quant] = obj.plotFavouriteStimulus(events, stim_types, cfg, nShuffles, binSize)
        %   [hf, tuningCurves, quant] = obj.plotFavouriteStimulus(events, stim_types, cfg, nShuffles, binSize, false)
        %
        %   INPUTS
        %     events      [T x N x trials]  peri-event activity (good neurons only)
        %     stim_types  [trials x 1]      trial stimulus labels (cell str or numeric)
        %     cfg         PlotConfig        optional, default PlotConfig()
        %     nShuffles   scalar            shuffles for null distribution (default 500)
        %     binSize     scalar            distance bin width in centroid units — px by
        %                                  default, µm when anatomy-aligned (default 20)
        %     doPlot      logical           produce figure (default true); false = quant only
        %
        %   OUTPUTS
        %     hf            figure handle (nReps+1 tiles)
        %     tuningCurves  [N x nStims x nReps]  per-cell per-rep tuning curves
        %     quant         struct with fields:
        %       .fav             [N x nReps]      favourite stimulus index per rep (NaN = invalid)
        %       .fav_avg         [N x 1]          favourite from rep-averaged tuning curve
        %       .rep_curves      [nBins x nReps]  P(same fav) binned by distance, per rep
        %       .avg_curve       [nBins x 1]      P(same fav) from rep-averaged map
        %       .rep_shuf_mean   [nBins x nReps]  per-rep shuffle mean
        %       .rep_shuf_std    [nBins x nReps]  per-rep shuffle std
        %       .avg_shuf_mean   [nBins x 1]      avg-map shuffle mean
        %       .avg_shuf_std    [nBins x 1]      avg-map shuffle std
        %       .binCenters      [1 x nBins]      bin centre distances
        %       .binEdges        [1 x nBins+1]
        %       .binSize         scalar            bin width (same units as centroids)
        %       .stim_labels     {nStims x 1}
            arguments
                obj
                events (:,:,:) double
                stim_types
                cfg PlotConfig = PlotConfig()
                nShuffles (1,1) double = 500
                binSize (1,1) double = 20
                doPlot (1,1) logical = true
            end

            N = numel(obj.roi_ids);
            assert(size(events, 2) == N, ...
                'events must have %d units (columns), got %d', N, size(events, 2));

            % ---- 1. Tuning curves: [N x nStims x nReps] ----
            tuningCurves = extractActivityMetric(events, 'tuning curves', 'cells', ...
                'StimTypes', stim_types);
            [~, nStims, nReps] = size(tuningCurves);

            if iscell(stim_types)
                stim_labels = unique(stim_types(:));
            else
                stim_labels = arrayfun(@num2str, unique(stim_types(:)), 'UniformOutput', false);
            end

            % ---- 2. Favourite stimulus per cell per rep (argmax over positive responses) ----
            fav = nan(N, nReps);
            for r = 1:nReps
                tc_r = tuningCurves(:, :, r);   % [N x nStims]
                [mx, idx] = max(tc_r, [], 2);
                valid = mx > 0;
                fav(valid, r) = idx(valid);
            end

            % ---- 3. Rep-averaged favourite ----
            tc_avg = mean(tuningCurves, 3, 'omitmissing');   % [N x nStims]
            [mx_avg, idx_avg] = max(tc_avg, [], 2);
            fav_avg = nan(N, 1);
            fav_avg(mx_avg > 0) = idx_avg(mx_avg > 0);

            % ---- 4. Pairwise anatomical distances (upper triangle) ----
            distMat  = squareform(pdist(obj.centroids, 'euclidean'));
            tri_mask = triu(true(N), 1);
            dists    = distMat(tri_mask);

            % Adaptive bins: fixed width, edges snapped to multiples of binSize
            d_min      = floor(min(dists) / binSize) * binSize;
            d_max      = ceil(max(dists)  / binSize) * binSize;
            binEdges   = d_min : binSize : d_max;
            binCenters = binEdges(1:end-1) + binSize / 2;
            nBins      = numel(binCenters);
            binIdx     = discretize(dists, binEdges);

            % ---- 5. Binned P(same fav) for each rep ----
            rep_curves = nan(nBins, nReps);
            for r = 1:nReps
                rep_curves(:, r) = computeSameFavCurve(fav(:, r), tri_mask, binIdx, nBins);
            end

            % ---- 6. Rep-averaged curve ----
            avg_curve = computeSameFavCurve(fav_avg, tri_mask, binIdx, nBins);

            % ---- 7. Per-curve label shuffles ----
            % One independent shuffle per curve (per rep + avg) so each can be
            % normalized against its own null, removing odour-dominance offsets.
            rep_shuf_mean = nan(nBins, nReps);
            rep_shuf_std  = nan(nBins, nReps);
            for r = 1:nReps
                vi_r         = find(~isnan(fav(:, r)));
                shuf_r       = nan(nBins, nShuffles);
                for s = 1:nShuffles
                    fs           = fav(:, r);
                    fs(vi_r)     = fav(vi_r(randperm(numel(vi_r))), r);
                    shuf_r(:, s) = computeSameFavCurve(fs, tri_mask, binIdx, nBins);
                end
                rep_shuf_mean(:, r) = mean(shuf_r, 2, 'omitmissing');
                rep_shuf_std(:, r)  = std(shuf_r,  [], 2, 'omitmissing');
            end

            vi_avg        = find(~isnan(fav_avg));
            shuf_avg      = nan(nBins, nShuffles);
            for s = 1:nShuffles
                fs            = fav_avg;
                fs(vi_avg)    = fav_avg(vi_avg(randperm(numel(vi_avg))));
                shuf_avg(:,s) = computeSameFavCurve(fs, tri_mask, binIdx, nBins);
            end
            avg_shuf_mean = mean(shuf_avg, 2, 'omitmissing');
            avg_shuf_std  = std(shuf_avg,  [], 2, 'omitmissing');

            % ---- 8. Package output ----
            quant = struct( ...
                'fav',           fav, ...
                'fav_avg',       fav_avg, ...
                'rep_curves',    rep_curves, ...
                'avg_curve',     avg_curve, ...
                'rep_shuf_mean', rep_shuf_mean, ...
                'rep_shuf_std',  rep_shuf_std, ...
                'avg_shuf_mean', avg_shuf_mean, ...
                'avg_shuf_std',  avg_shuf_std, ...
                'binCenters',    binCenters, ...
                'binEdges',      binEdges, ...
                'binSize',       binSize, ...
                'stim_labels',   {stim_labels});

            % ---- 9. Plot (optional) ----
            hf = gobjects(0);
            if ~doPlot; return; end

            % Stimulus colormap (categorical, one color per stimulus)
            if nStims <= size(cfg.c, 1)
                stim_cmap = cfg.c(1:nStims, :);
            else
                stim_cmap = lines(nStims);
            end

            % Figure: nReps topo maps + 1 quantification tile
            hf = figure('Color', cfg.bgcol);
            tl = tiledlayout(hf, 1, nReps + 1, 'TileSpacing', 'compact', 'Padding', 'compact');

            for r = 1:nReps
                ax = nexttile(tl, r);
                obj.plotTopoMap(fav(:, r), ax);
                colormap(ax, stim_cmap);
                clim(ax, [0.5, nStims + 0.5]);
                cb = ax.Colorbar;
                cb.Ticks = 1:nStims;
                cb.TickLabels = stim_labels;
                cb.Label.String = 'Fav. stimulus';
                cb.Label.Color  = cfg.textcol;
                cb.Color        = cfg.axcol;
                title(ax, sprintf('Rep %d', r), 'Color', cfg.textcol);
            end

            % Quantification tile (nReps+2 distributions via shared utility)
            nexttile(tl, nReps + 1);
            plotFavStimulusCurves(quant, cfg, gca);
        end

        function results = quantify(obj, inMat, similarityMetric, nShuffles, nBins)
            arguments
                obj
                inMat double          % (T, units)
                similarityMetric char = 'correlation'
                nShuffles double = 1000
                nBins double = 10
            end
            N = numel(obj.roi_ids);
            assert(size(inMat, 2) == N, ...
                'inMat must have %d columns (one per good neuron)', N);

            % pairwise anatomical distance
            distVec = pdist(obj.centroids, 'euclidean');
            distMat = squareform(distVec);

            % pairwise functional similarity
            simVec = computeSimilarity(inMat, similarityMetric);
            simMat = squareform(simVec);

            % upper triangle indices
            mask = triu(true(N), 1);
            dists = distMat(mask);
            sims  = simMat(mask);

            % uniform distance bins
            binEdges = linspace(min(dists), max(dists) + eps, nBins + 1);
            binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2;
            [binIdx, ~] = discretize(dists, binEdges);

            % real binned curve
            realCurve = zeros(nBins, 1);
            realSEM   = zeros(nBins, 1);
            for b = 1:nBins
                vals = sims(binIdx == b);
                realCurve(b) = mean(vals, 'omitmissing');
                realSEM(b)   = std(vals, 'omitmissing') / sqrt(sum(~isnan(vals)));
            end

            % Mantel test statistic (real)
            mantelR_real = corr(dists, sims, 'rows', 'pairwise');

            % shuffle
            shuffCurves = zeros(nBins, nShuffles);
            mantelR_shuf = zeros(nShuffles, 1);
            for s = 1:nShuffles
                perm = randperm(N);
                simVec_s = computeSimilarity(inMat(:, perm), similarityMetric);
                simMat_s = squareform(simVec_s);
                sims_s = simMat_s(mask);
                for b = 1:nBins
                    shuffCurves(b, s) = mean(sims_s(binIdx == b), 'omitmissing');
                end
                mantelR_shuf(s) = corr(dists, sims_s, 'rows', 'pairwise');
            end

            % shuffle statistics
            shuffMean = mean(shuffCurves, 2);
            shuffCI   = prctile(shuffCurves, [2.5 97.5], 2); % (nBins, 2)

            % per-bin p-values (two-tailed)
            pBin = zeros(nBins, 1);
            for b = 1:nBins
                pBin(b) = mean(abs(shuffCurves(b,:)) >= abs(realCurve(b)));
            end

            % Mantel p-value (two-tailed)
            mantelP = mean(abs(mantelR_shuf) >= abs(mantelR_real));

            % plot
            hf = figure; hold on;
            cfg.lineStyle = '-'; cfg.lineWidth = 1;

            % shuffle ribbon
            plotLineNShade(binCenters, shuffMean, ...
                (shuffCI(:,2) - shuffCI(:,1)) / 2, [0.6 0.6 0.6], cfg);

            % real data
            cfg.lineWidth = 2;
            plotLineNShade(binCenters, realCurve, realSEM, [0 0 0], cfg);

            xlabel('Anatomical distance (px)');
            ylabel(['Functional similarity (' similarityMetric ')']);
            title(sprintf('Mantel r = %.3f, p = %.4f', mantelR_real, mantelP));
            legend({'shuffle 95% CI', '', 'data', ''}, 'Location', 'best');
            hold off;

            % output
            results.binCenters = binCenters;
            results.binEdges   = binEdges;
            results.realCurve  = realCurve;
            results.realSEM    = realSEM;
            results.shuffMean  = shuffMean;
            results.shuffCI    = shuffCI;
            results.pBin       = pBin;
            results.mantelR    = mantelR_real;
            results.mantelP    = mantelP;
            results.hf         = hf;
        end
    end
end

function curve = computeSameFavCurve(fav_vec, tri_mask, binIdx, nBins)
% Per-distance-bin fraction of cell pairs sharing the same favourite stimulus.
% fav_vec: [N x 1] double, NaN where no valid favourite was assigned.
    N    = numel(fav_vec);
    vi   = find(~isnan(fav_vec));

    same_mat = nan(N, N);
    if numel(vi) >= 2
        sub             = fav_vec(vi) == fav_vec(vi)';   % [nValid x nValid]
        same_mat(vi,vi) = double(sub);
    end

    same_vec = same_mat(tri_mask);
    curve    = nan(nBins, 1);
    for b = 1:nBins
        vals = same_vec(binIdx == b);
        vals = vals(~isnan(vals));
        if ~isempty(vals)
            curve(b) = mean(vals);
        end
    end
end

function simVec = computeSimilarity(inMat, metric)
    valid = sum(ismissing(inMat),2)==0;
    d = pdist(inMat(valid,:)', metric);
    switch metric
        case 'correlation'
            simVec = 1 - d; % correlation distance → similarity
        case 'cosine'
            simVec = 1 - d;
        otherwise
            simVec = -d; % for euclidean etc., negate so higher = more similar
    end
end
