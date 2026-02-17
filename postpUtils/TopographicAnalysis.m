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
