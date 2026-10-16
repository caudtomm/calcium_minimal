function [b, params] = prettyBoxplot(dataMat, labs, varargin)
arguments
    dataMat double
    labs cell = {''}
end
arguments (Repeating)
    varargin
end

params = paramsStruct();
for k = 1:2:length(varargin)
    paramName = varargin{k};
    paramValue = varargin{k+1};
    if isfield(params, paramName)
        params.(paramName) = paramValue;
    else
        error('Unknown parameter name: %s', paramName);
    end

end

[nvars, nmeas] = size(dataMat);

hold on
datacells = num2cell(dataMat, 1); % {1 x nmeas}, each [nvars x 1]
% box = IQR, median line, optional mean+/-SD whiskers, optional mean marker
RF_mkBoxPlot3(datacells, params.boxColor, 1:nmeas, params.boxWidth, ...
    params.showSD, params.boxLineWidth, params.meanMarkerSize, []);

xticks(1:nmeas);
xticklabels(labs);
xlim([0.5, nmeas+0.5]);

b(1) = scatter(1:nmeas,mean(dataMat,'omitmissing'), ...
    params.scatterSize, ...
    params.scatterColor,'filled');

box off;

if ~params.plotLine; return; end

b(2) = plot([0 nmeas+1],params.lineVal * [1 1], ...
    'LineWidth', params.lineWidth, ...
    'Color', params.lineColor, ...
    'LineStyle', params.lineStyle);

end

function params = paramsStruct()
params = struct(...
    'boxColor', [0 0 0], ...
    'boxWidth', [], ...       % [] -> RF_mkBoxPlot3 default width (0.8)
    'boxLineWidth', 1, ...
    'showSD', true, ...      % overlay mean +/- 1 SD whiskers, RF_mkBoxPlot3-style
    'meanMarkerSize', 0, ...  % RF_mkBoxPlot3's own mean marker (box-colored); left off by default, see scatter* below
    ...
    'scatterColor', [1 0 0], ...
    'scatterSize', 36, ...
    'scatterMarker', 'o', ...
    ...
    'plotLine', false, ...
    'lineVal', 0, ...
    'lineWidth', 1.5, ...
    'lineColor', [1 0 0], ...
    'lineStyle', '--' ...
    );
end