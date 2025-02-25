function plotHeatmapAndIsoclines(x, y, gridSize)
    % % Create a grid for the heatmap
    xEdges = linspace(min(x), max(x), gridSize);
    yEdges = linspace(min(y), max(y), gridSize);
    
    % Compute a 2D histogram (density map)
    density = histcounts2(x, y, xEdges, yEdges);
    
    % Convert edges to bin centers
    xCenters = xEdges(1:end-1) + diff(xEdges)/2;
    yCenters = yEdges(1:end-1) + diff(yEdges)/2;
    
    % Create the heatmap
    figure;
    imagesc(xCenters, yCenters, density.');
    axis xy;
    colorbar;
    title('Heatmap of Point Density');
    xlabel('X');
    ylabel('Y');
    
    % Overlay isoclines (contours) on the heatmap
    hold on;
    [X, Y] = meshgrid(xCenters, yCenters);
    contour(X, Y, density.', 'LineColor', 'k', 'LineWidth', 1.5);
    hold off;
end