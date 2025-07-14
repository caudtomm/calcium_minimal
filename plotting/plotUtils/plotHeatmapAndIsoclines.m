function plotHeatmapAndIsoclines(x, y, gridSize, do_heatmap, do_isoclines, do_log)
    % % Create a grid for the heatmap
    xEdges = linspace(min(x), max(x), gridSize);
    yEdges = linspace(min(y), max(y), gridSize);
    
    % Compute a 2D histogram (density map)
    if ~do_log
        density = histcounts2(x, y, xEdges, yEdges);
    else
        density = log(histcounts2(x, y, xEdges, yEdges));
    end
    
    % Convert edges to bin centers
    xCenters = xEdges(1:end-1) + diff(xEdges)/2;
    yCenters = yEdges(1:end-1) + diff(yEdges)/2;
    
    % Create the heatmap
    if do_heatmap
        imagesc(xCenters, yCenters, density.');
        axis xy;
        colorbar;
        title('Heatmap of Point Density');
    end
    
    % Overlay isoclines (contours) on the heatmap
    if do_isoclines
        hold on;
        [X, Y] = meshgrid(xCenters, yCenters);
        contour(X, Y, density.', 'LineColor', 'k', 'LineWidth', 1.5);
        hold off;
    end

    xlabel('X');
    ylabel('Y');
end