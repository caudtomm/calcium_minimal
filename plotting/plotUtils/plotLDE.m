function out = plotLDE(embedding,plotType,labs,cfg, varargin)
% <Low level low dimensional embedding plotter>
arguments
    embedding double % [nSamples x nDimensions x nTrials] matrix
    plotType char = 'scatter'
    labs = [] % stimulus names (sorted!)
    cfg = PlotConfig()
end
arguments (Repeating)
    varargin
end

% Default parameters
params.ldeType = 'projection';
params.sizeByRep = true; % size points by repetition number
params.showNumbers  = true;  % trialNum: plot absolute trial number at COM
params.showTrialLine = true; % trialNum: blue line connecting consecutive trial COMs
params.showRepLines  = true; % trialNum: red lines connecting same-stimulus COMs

% Parse name-value pairs
if ~isempty(varargin)
    for k = 1:2:length(varargin)
        switch lower(varargin{k})
            case 'ldetype'
                params.ldeType = varargin{k+1};
            case 'shownumbers'
                params.showNumbers = varargin{k+1};
            case 'showtrialline'
                params.showTrialLine = varargin{k+1};
            case 'showreplines'
                params.showRepLines = varargin{k+1};
            otherwise
                error('Unknown parameter name: %s', varargin{k});
        end
    end
end

% Parse input
assert(numel(labs) == size(embedding, 3), 'Mismatch between labs and trial dimension.'); % labs length validation

% Initialize output
out=[];

%% Plot onto provided axes
switch plotType
    case 'scatter' % Plot all in one big matrix
        out = plotScatter(embedding,labs,params,cfg);
    case 'lines'
        out = plotLines(embedding,labs,params,cfg);
    case 'trialNum'
        out = plotTrialNum(embedding,labs,params,cfg);
    otherwise
        error('Requested plot type is unknown.')
end

end


%% Functions
function h = plotScatter(embedding,labs,params,cfg)
    [L, ndims, ntrials] = size(embedding);

    % set scatter size
    scalefact = 500;
    sizesc = 10+scalefact/(L*ntrials);

    stims = unique(labs,'stable');

    h = gobjects(1, ntrials); % Preallocate h array for better performance
    for i_trial = 1:ntrials
        thisStim = labs(i_trial);
        thisStimID = find(ismember(stims,thisStim));
        thisRepNum = sum(ismember(labs(1:i_trial),thisStim));

        if numel(stims) > 1
            thisColor = cfg.c(thisStimID,:); % color by stimulus
        else
            thisColor = cfg.c(thisRepNum,:); %color by repetition
        end

        sizeMultiplier = 1;
        if params.sizeByRep && numel(stims) > 1
            sizeMultiplier = 1/thisRepNum;
        end

        switch ndims
            case 2
                h(i_trial) = scatter(embedding(:,1,i_trial),embedding(:,2,i_trial), ...
                    sizesc*sizeMultiplier,thisColor,'filled');
            case 3
                h(i_trial) = scatter3(embedding(:,1,i_trial),embedding(:,2,i_trial), ...
                    embedding(:,3,i_trial),sizesc*sizeMultiplier, ...
                    thisColor,'filled');
            otherwise
                error('dimension number not supported')
        end

        if numel(stims) > 1
            h(i_trial).DisplayName = thisStim{1}; % legend by stimulus
        else
            h(i_trial).DisplayName = ['Rep #',num2str(thisRepNum)]; % legend by repetition
        end
        hold on
    end


    % cosmetics / labels
    u = legendUnq();
    legend(u,'Box','on','color',cfg.bgcol,'Location','best', ...
        'EdgeColor',cfg.textcol,'TextColor',cfg.textcol)
    set(gcf, 'color', cfg.bgcol);
    set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, ...
        'YColor',cfg.axcol, 'ZColor',cfg.axcol);
    xlabel([params.ldeType,' #1'])
    ylabel([params.ldeType,' #2'])
    if ndims>2; zlabel([params.ldeType,' #3']); end
    grid off
    axis square
end

function h = plotLines(embedding,labs,params,cfg)
    [L, ndims, ntrials] = size(embedding);

    % set scatter size
    scalefact = 5;
    sizesc = 10+scalefact/(L*ntrials);

    stims = unique(labs, 'stable');

    for i_trial = 1:ntrials
        thisStim = labs(i_trial);
        thisStimID = find(ismember(stims,thisStim));
        thisRepNum = sum(ismember(labs(1:i_trial),thisStim));

        thisColor = cfg.c(thisStimID,:);

        sizeMultiplier = 1;
        if params.sizeByRep
            sizeMultiplier = .3/thisRepNum;
        end

        switch ndims
            case 2
                h(i_trial) = plot(embedding(:,1,i_trial),embedding(:,2,i_trial), ...
                    'LineWidth',sizesc*sizeMultiplier,'Color',thisColor);
            case 3
                h(i_trial) = plot3(embedding(:,1,i_trial),embedding(:,2,i_trial), ...
                    embedding(:,3,i_trial), ...
                    'LineWidth',sizesc*sizeMultiplier, ...
                    'Color',thisColor);
            otherwise
                error('dimension number not supported')
        end
        h(i_trial).DisplayName = thisStim{1};
        hold on
    end


    % cosmetics / labels
    u = legendUnq();
    legend(u,'Box','on','color',cfg.bgcol,'Location','best', ...
        'EdgeColor',cfg.textcol,'TextColor',cfg.textcol)
    set(gcf, 'color', cfg.bgcol);
    set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, ...
        'YColor',cfg.axcol, 'ZColor',cfg.axcol);
    xlabel([params.ldeType,' #1'])
    ylabel([params.ldeType,' #2'])
    if ndims>2; zlabel([params.ldeType,' #3']); end
    grid off
    axis square
end

function h = plotTrialNum(embedding,labs,params,cfg)
    [~, ndims, ntrials] = size(embedding);
    stims = unique(labs, 'stable');

    if ~params.showNumbers && ~params.showTrialLine && ~params.showRepLines
        error('plotLDE:trialNum:nothingToPlot', ...
            'At least one of showNumbers, showTrialLine, or showRepLines must be true.');
    end

    % Compute center of mass for every trial upfront
    coms = zeros(ntrials, ndims);
    for i = 1:ntrials
        coms(i,:) = extractCenterOfMass(embedding, i);
    end

    h = struct('numbers', {{}}, 'trial_line', {{}}, 'rep_lines', {{}});

    % Linewidth decreases linearly from max to min with increasing trial number
    lw_scale = @(trial_num) max(0.2, 1.5 - 1.2 * (trial_num - 1) / max(ntrials - 1, 1));

    % Blue semitransparent line connecting consecutive trial COMs
    if params.showTrialLine && ntrials > 1
        for i = 1:ntrials-1
            hl = plotConnectingLine(coms(i,:), coms(i+1,:), ndims, ...
                'Color', [0.2 0.4 1 0.35], 'LineWidth', lw_scale(i));
            hold on
            h.trial_line{end+1} = hl;
        end
    end

    % Red lines per stimulus connecting consecutive same-stimulus COMs
    if params.showRepLines
        for i_stim = 1:numel(stims)
            stim_idx = find(ismember(labs, stims(i_stim)));
            for k = 1:numel(stim_idx)-1
                hl = plotConnectingLine(coms(stim_idx(k),:), coms(stim_idx(k+1),:), ndims, ...
                    'Color', [1 0.2 0.2 0.5], 'LineWidth', lw_scale(stim_idx(k)) * 1.2);
                hold on
                h.rep_lines{end+1} = hl;
            end
        end
    end

    % Absolute trial number text at each COM, colored by stimulus
    if params.showNumbers
        for i_trial = 1:ntrials
            thisStim = labs(i_trial);
            thisStimID = find(ismember(stims, thisStim));
            if numel(stims) > 1
                thisColor = cfg.c(thisStimID,:);
            else
                thisRepNum = sum(ismember(labs(1:i_trial), thisStim));
                thisColor = cfg.c(thisRepNum,:);
            end
            ht = plotNumberAtPoint(coms(i_trial,:), i_trial, ndims, thisColor);
            hold on
            h.numbers{end+1} = ht;
        end
    end

    % cosmetics / labels
    set(gcf, 'color', cfg.bgcol);
    set(gca, 'color', cfg.bgcol, 'XColor', cfg.axcol, ...
        'YColor', cfg.axcol, 'ZColor', cfg.axcol);
    xlabel([params.ldeType,' #1'])
    ylabel([params.ldeType,' #2'])
    if ndims > 2; zlabel([params.ldeType,' #3']); end
    grid off
    axis square
end

% --- shared low-level helpers ---

function com = extractCenterOfMass(embedding, i_trial)
    com = mean(embedding(:,:,i_trial), 1); % [1 x ndims]
end

function h = plotConnectingLine(p1, p2, ndims, varargin)
    switch ndims
        case 2
            h = plot([p1(1) p2(1)], [p1(2) p2(2)], varargin{:});
        case 3
            h = plot3([p1(1) p2(1)], [p1(2) p2(2)], [p1(3) p2(3)], varargin{:});
        otherwise
            error('dimension number not supported')
    end
end

function h = plotNumberAtPoint(pt, num, ndims, color)
    label = num2str(num);
    opts = {'Color', color, 'FontSize', 9, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle'};
    switch ndims
        case 2
            h = text(pt(1), pt(2), label, opts{:});
        case 3
            h = text(pt(1), pt(2), pt(3), label, opts{:});
        otherwise
            error('dimension number not supported')
    end
end
