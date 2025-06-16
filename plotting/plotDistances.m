function plotDistances(data)
arguments
    data struct
end

% --- Defaults ---
defaults = struct( ...
    'colormap', parula, ...
    'axisColor', [0 0 0], ...
    'textColor', [0 0 0], ...
    'backgroundColor', [1 1 1], ...
    'lineWidth', 1.5 ...
);

cfg = parseThemeColors(data,defaults,'themeColors');




end