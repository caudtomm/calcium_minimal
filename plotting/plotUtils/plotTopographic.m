function out = plotTopographic(all_metrics, all_traces, plotType, events, cfg, varargin)
% PLOTTOPOGRAPHIC Low-level topographic analysis plotter.
%
%   out = plotTopographic(all_metrics, all_traces, plotType, events, cfg)
%
%   INPUTS:
%       all_metrics  - cell array {nSubjects x 1}, each cell is [units x 1] metric
%       all_traces   - cell array {nSubjects x 1} of ActivityTraces objects
%       plotType     - 'maps' or 'stats'
%       events       - cell array {nSubjects x 1}, each [T x units x trials]
%                      (used only for 'stats' mode)
%       cfg          - PlotConfig object
%
%   OPTIONAL NAME-VALUE PAIRS:
%       'similarityMetric' - metric for quantify() (default: 'correlation')
%       'nShuffles'        - number of shuffles (default: 1000)
%       'nBins'            - number of distance bins (default: 10)
%
%   OUTPUTS:
%       out - struct with figure handles and results

% parse optional args
similarityMetric = 'correlation';
nShuffles = 1000;
nBins = 10;
for k = 1:2:numel(varargin)
    switch lower(varargin{k})
        case 'similaritymetric'
            similarityMetric = varargin{k+1};
        case 'nshuffles'
            nShuffles = varargin{k+1};
        case 'nbins'
            nBins = varargin{k+1};
    end
end

nsubjects = numel(all_metrics);
out = struct();

% build TopographicAnalysis objects
tas = cell(nsubjects, 1);
for i = 1:nsubjects
    tas{i} = TopographicAnalysis(all_traces{i});
end

% normalize metrics per subject: z-score
norm_metrics = cell(nsubjects, 1);
for i = 1:nsubjects
    m = all_metrics{i};
    norm_metrics{i} = (m - mean(m, 'omitmissing')) / std(m, 'omitmissing');
end

switch lower(plotType)
    case 'maps'
        % --- Figure 1: individual subject topo maps ---
        ncols = ceil(sqrt(nsubjects));
        nrows = ceil(nsubjects / ncols);
        hf1 = figure;
        for i = 1:nsubjects
            ax = subplot(nrows, ncols, i);
            tas{i}.plotTopoMap(norm_metrics{i}, ax);
            title(ax, sprintf('Subject %d', i), 'Color', cfg.textcol);
        end
        set(hf1, 'Color', cfg.bgcol);
        out.hf_individual = hf1;

        % --- Figure 2: overlay scatter of all centroids ---
        hf2 = figure; hold on;
        all_x = []; all_y = []; all_v = [];
        for i = 1:nsubjects
            c = tas{i}.centroids;
            all_x = [all_x; c(:,1)]; %#ok<AGROW>
            all_y = [all_y; c(:,2)]; %#ok<AGROW>
            all_v = [all_v; norm_metrics{i}]; %#ok<AGROW>
        end
        scatter(all_x, all_y, 30, all_v, 'filled', 'CData', all_v);
        colormap(gca, 'parula'); axis square; clim([-1 1])
        colorbar;
        set(gca, 'YDir', 'reverse', 'Color', cfg.bgcol, ...
            'XColor', cfg.axcol, 'YColor', cfg.axcol);
        axis equal tight; box off;
        xlabel('x (px)'); ylabel('y (px)');
        title('All subjects (normalized)', 'Color', cfg.textcol);
        set(hf2, 'Color', cfg.bgcol);
        hold off;
        out.hf_overlay = hf2;

        
        % --- Figure 3: spatial heatmap of interpolated metric values ---
        hf3 = figure;
        
        % spatial grid
        gridSize = 50;
        xEdges = linspace(min(all_x), max(all_x), gridSize + 1);
        yEdges = linspace(min(all_y), max(all_y), gridSize + 1);
        [XI, YI] = meshgrid((xEdges(1:end-1) + xEdges(2:end)) / 2, ...
                    (yEdges(1:end-1) + yEdges(2:end)) / 2);
        
        % Interpolate metric values onto grid using scatteredInterpolant
        F = scatteredInterpolant(all_x, all_y, all_v, 'linear', 'none');
        ZI = F(XI, YI);
        
        imagesc(xEdges([1 end]), fliplr(yEdges([1 end])), ZI);
        axis square; clim([-1 1])
        set(gca, 'YDir', 'normal');
        colormap(gca, 'parula');
        colorbar;
        set(gca, 'Color', cfg.bgcol, 'XColor', cfg.axcol, 'YColor', cfg.axcol);
        xlabel('x (px)'); ylabel('y (px)');
        title('Spatial heatmap (normalized)', 'Color', cfg.textcol);
        set(hf3, 'Color', cfg.bgcol);
        out.hf_heatmap = hf3;
    case 'stats' % the pooled quantification doesn't make sense. the cells are not actually close to each other if they are from different subjects
        % pool activity data and centroids across subjects, after normalizing
        all_inMat = [];
        all_centroids = [];
        for i = 1:nsubjects
            inMat_i = events{i}; % [T x units x trials]
            if ndims(inMat_i) == 3
                inMat_i = ActivityTraces.format(inMat_i);
            end
            % z-score each unit's timecourse within subject
            inMat_i = nanzscore(inMat_i);
            all_inMat = [all_inMat, inMat_i]; %#ok<AGROW>
            all_centroids = [all_centroids; tas{i}.centroids]; %#ok<AGROW>
        end

        % build a temporary TopographicAnalysis with pooled data
        pooled = struct();
        pooled.centroids = all_centroids;
        pooled.roi_ids = 1:size(all_centroids, 1); % dummy ids

        % call quantify logic directly (avoid needing a full ActivityTraces)
        results = quantifyPooled(all_inMat, all_centroids, ...
            similarityMetric, nShuffles, nBins, cfg);
        out.results = results;

    otherwise
        error('Unknown plotType: %s. Use ''maps'' or ''stats''.', plotType);
end

end


function results = quantifyPooled(inMat, centroids, similarityMetric, nShuffles, nBins, cfg)
    N = size(inMat, 2);

    % pairwise anatomical distance
    distVec = pdist(centroids, 'euclidean');
    dists = squareform(distVec);

    % pairwise functional similarity
    simVec = computeSimVec(inMat, similarityMetric);
    sims = squareform(simVec);

    % upper triangle
    mask = triu(true(N), 1);
    d = dists(mask);
    s = sims(mask);

    % uniform distance bins
    binEdges = linspace(min(d), max(d) + eps, nBins + 1);
    binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2;
    binIdx = discretize(d, binEdges);

    % real binned curve
    realCurve = zeros(nBins, 1);
    realSEM   = zeros(nBins, 1);
    for b = 1:nBins
        vals = s(binIdx == b);
        realCurve(b) = mean(vals, 'omitmissing');
        realSEM(b)   = std(vals, 'omitmissing') / sqrt(sum(~isnan(vals)));
    end

    % Mantel test
    mantelR_real = corr(d, s, 'rows', 'pairwise');

    % shuffle
    shuffCurves = zeros(nBins, nShuffles);
    mantelR_shuf = zeros(nShuffles, 1);
    for sh = 1:nShuffles
        perm = randperm(N);
        simVec_s = computeSimVec(inMat(:, perm), similarityMetric);
        sims_s = squareform(simVec_s);
        s_s = sims_s(mask);
        for b = 1:nBins
            shuffCurves(b, sh) = mean(s_s(binIdx == b), 'omitmissing');
        end
        mantelR_shuf(sh) = corr(d, s_s, 'rows', 'pairwise');
    end

    shuffMean = mean(shuffCurves, 2);
    shuffCI   = prctile(shuffCurves, [2.5 97.5], 2);

    pBin = zeros(nBins, 1);
    for b = 1:nBins
        pBin(b) = mean(abs(shuffCurves(b,:)) >= abs(realCurve(b)));
    end
    mantelP = mean(abs(mantelR_shuf) >= abs(mantelR_real));

    % plot
    hf = figure; hold on;
    lcfg.lineStyle = '-'; lcfg.lineWidth = 1;
    plotLineNShade(binCenters, shuffMean, ...
        (shuffCI(:,2) - shuffCI(:,1)) / 2, [0.6 0.6 0.6], lcfg);
    lcfg.lineWidth = 2;
    plotLineNShade(binCenters, realCurve, realSEM, cfg.c(1,:), lcfg);

    xlabel('Anatomical distance (px)');
    ylabel(['Functional similarity (' similarityMetric ')']);
    title(sprintf('Pooled: Mantel r = %.3f, p = %.4f', mantelR_real, mantelP), ...
        'Color', cfg.textcol);
    legend({'shuffle 95% CI', '', 'data', ''}, 'Location', 'best');
    set(gca, 'Color', cfg.bgcol, 'XColor', cfg.axcol, 'YColor', cfg.axcol);
    set(hf, 'Color', cfg.bgcol);
    hold off;

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


function simVec = computeSimVec(inMat, metric)
    d = pdist(inMat', metric);
    switch metric
        case 'correlation'
            simVec = 1 - d;
        case 'cosine'
            simVec = 1 - d;
        otherwise
            simVec = -d;
    end
end
