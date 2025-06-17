function [hf, out] = plotDistances(data)
arguments
    data struct
end

% Knobs
method = 'correlation';
crange = [.2 1];

% Input data parsing
% expected input:
% .traces - cell array [nsubjects 1] of double [t, cells, trials] (sorted!)
% .fs - framerate in seconds (assumes the same for all)
% .sec_range - double [absolute absolute] (assumes the same for all)
% .labs - cell array [nsubjects 1] of double indices [ntrials 1] with
%       stimulus names (sorted!)
traces = data.traces;
fs = data.fs;
sec_range = data.sec_range;
labs = data.labs;


% Color Theme Configuration parsing

% --- Defaults ---
defaults = struct( ...
    'colormap', parula, ...
    'axisColor', [0 0 0], ...
    'textColor', [0 0 0], ...
    'backgroundColor', [1 1 1], ...
    'lineWidth', 1.5 ...
);

cfg = parseThemeColors(data,defaults,'themeColors');


% initialize useful metrics
nsubjects = numel(traces);
ntrials = numel(labs);


% Initialize output
hf = gobjects(1,1);
out.distanceMAT=[];


% Plot all in one big matrix
[hf(1), out.distanceMAT] = plotFullMat(true);






function [hf, distanceMAT3D] = plotFullMat(pl)
    distanceMAT3D = nan(ntrials,ntrials,nsubjects);
    interval = [floor(sec_range(1)*fs) : ...
              1+floor(sec_range(2)*fs)];
    for i_fish = 1:nsubjects
        c = calcdistanceMAT(traces{i_fish},interval,method);
        
        % store to common matrix
        distanceMAT3D(1:size(c,1),1:size(c,2),i_fish) = c; % indices are here in case this fish has less than a full trial set
    end

    % (optional) plot full avg intertrial distance matrix
    if ~pl; return; end
    distanceMAT2D = mean(distanceMAT3D,3,'omitmissing');
    hf = figure; imagesc(1-distanceMAT2D)
    axis square; hold on
    xticks(1:ntrials); xticklabels(labs); xtickangle(90)
    yticks(1:ntrials); yticklabels(labs)
    xlabel('Stimulus type')
    ylabel('Stimulus type')
    title(['Similarity: ',method, sec_range(1),'-',sec_range(2), ' s'], 'Color',cfg.textColor)
    clim(crange)
    colormap(cfg.colormap)
    colorbar('Color',cfg.axisColor)
    set(gca, 'color', cfg.backgroundColor, 'XColor',cfg.axisColor, 'YColor',cfg.axisColor, 'ZColor',cfg.axisColor);
    set(gcf, 'color', cfg.backgroundColor); 
    hold off
end

end

function distanceMAT = calcdistanceMAT(traces,interval,method)
snippets = traces(interval,:,:);
avg_actvect = squeeze(mean(snippets,1,'omitmissing'));
distanceMAT = squareform(pdist(avg_actvect', method));
end