function [vals, varargout] = extractActivityMetric(events, metric, n_equals, varargin)
% This function extracts various activity metrics from neural event data.
% 
% INPUTS:
%   events      - A 3D matrix of neural activity data [time x units x events/trials].
%   metric      - A string specifying the metric to calculate (e.g., 'population sparseness', 'tuning curves').
%   n_equals    - A string specifying the output format ('cells' or 'frames').
%   varargin    - Optional name-value pair arguments:
%                 'StimTypes': A vector specifying stimulus types for each trial.
%
% OUTPUTS:
%   vals        - The computed metric values. The size and format depend on the selected metric and 'n_equals'.
%   varargout   - Additional outputs, if applicable (e.g., average tuning curves for 'tuning curves' metric).
%
% DESCRIPTION:
%   This function processes neural activity data to compute various metrics
%   such as sparseness, tuning curves, selectivity, suppression scores, and more.
%   The input data is expected to be a 3D matrix where dimensions represent
%   time, units (neurons), and events/trials. The output format and size
%   depend on the selected metric and the 'n_equals' argument.
%
% EXAMPLES:
%   [vals] = extractActivityMetric(events, 'population sparseness', 'cells');
%   [vals] = extractActivityMetric(events, 'tuning curves', 'cells', 'StimTypes', stim_types);
%
% NOTES:
%   - The function includes several helper functions for specific calculations.
%   - Ensure that the 'metric' and 'n_equals' arguments are compatible.
%   - Some metrics require additional inputs (e.g., 'StimTypes').

% Validate inputs
if nargin < 3
    error('Not enough input arguments.');
end

% Parse name-value pair inputs
stim_types = [];
if ~isempty(varargin)
    for i = 1:2:numel(varargin)
        if strcmpi(varargin{i}, 'StimTypes')
            stim_types = varargin{i+1};
        end
    end
end

nSubjects = numel(events);
vals = [];

    data = events; % [time x units x events/trials]
    [L, nUnits, nEvents] = size(data);
    switch lower(n_equals)
        case 'cells'
            % Output: [units x events/trials]
            % Input : [time x units x events/trials] <-- nothing to do here
        case 'frames'
            % Output: [time x events/trials]
            % Input : [units x time x events/trials]
            data = permute(data, [2, 1, 3]); % [units x time x events/trials]
        otherwise
            error('Unknown n_equals: %s', n_equals);
    end
    
    switch lower(metric)

        % metrics that return a column vector of size [trials x 1]
        case 'population sparseness'
            checkn_equals(n_equals, 'cells');
            data = squeeze(mean(data,1,'omitmissing')); % get mean activity per unit across trials
            thisvals = squeeze(sum(data.^2, 1, 'omitmissing')) ./ ...
                (nUnits * sum(data, 1, 'omitmissing').^2);
        case 'normalized population sparseness'
            checkn_equals(n_equals, 'cells');
            data = squeeze(mean(data,1,'omitmissing')); % get mean activity per unit across trials
            thisvals = (squeeze(sum(data.^2, 1, 'omitmissing')) ./ ...
                (sum(data, 1, 'omitmissing').^2) - 1/nUnits) ./ (1 - 1/nUnits);
        case 'participation ratio'
            thisvals = zeros(1, nEvents);
            for i = 1:nEvents
                c = cov(data(:,:,i), 'omitrows'); % covariance matrix
                eigv = eig(c); % eigenvalues
                rel_eigv = eigv / sum(eigv); % normalized eigenvalues
                thisvals(i) = 1 ./ sum(rel_eigv.^2); % participation ratio
            end
       
        % metrics that return a column vector of size [units x stimuli]
        case 'tuning curves'
            checkn_equals(n_equals, 'cells');
            thisvals = getTuningCurves(data, stim_types); % [units x stimuli x repetitions]


        % metrics that return a column vector of size [units x repetitions]
        case 'selectivity of tuning repetitions'
            checkn_equals(n_equals, 'cells');
            checkexists(stim_types, 'StimTypes');

            % Calculate mean activity level for each neuron for each stimulus
            meanActivity = getTuningCurves(data, stim_types);

            % Calculate and store tuning selectivity for each neuron
            thisvals = calculateTuningSelectivity(meanActivity);


        % metrics that return a column vector of size [units x 1]
        case 'selectivity of tuning'
            checkn_equals(n_equals, 'cells');
            checkexists(stim_types, 'StimTypes');

            % Calculate mean activity level for each neuron for each stimulus
            [~, meanActivity] = getTuningCurves(data, stim_types);

            % Calculate and store tuning selectivity for each neuron
            thisvals = calculateTuningSelectivity(meanActivity);

        case 'stimulus specific suppression score'
            % answers the question: how much is the activity of a neuron 
            % suppressed over repetitions for one stimulus over all others?
            checkn_equals(n_equals, 'cells');
            checkexists(stim_types, 'StimTypes');
            
            % Calculate mean activity level for each neuron for each stimulus
            meanActivity = getTuningCurves(data, stim_types); % [numNeurons x numStims x numReps]

            suppression = getSuppressionScores(meanActivity); % [numNeurons x numStims]
            [numNeurons, nStims] = size(suppression);
            
            % compare across stimuli: how many units of STD of the others is the maximum away?
            suppression = sort(suppression, 2, 'descend'); % sort by suppression score
            suppression = suppression - mean(suppression(2:nStims), 2, 'omitnan'); % subtract mean suppression score of the non-maximum stimuli
            
            thisvals = suppression(:,1) / std(suppression(:,2:nStims), 0, 2, 'omitnan'); % divide by STD of the non-maximum stimuli

        case 'general suppression score'
            % answers the question: how much is the activity of a neuron 
            % suppressed over repetitions on average?
            checkn_equals(n_equals, 'cells');
            checkexists(stim_types, 'StimTypes');
            
            % Calculate mean activity level for each neuron for each stimulus
            meanActivity = getTuningCurves(data, stim_types); % [numNeurons x numStims x numReps]

            suppression = getSuppressionScores(meanActivity); % [numNeurons x numStims]
            [numNeurons, nStims] = size(suppression);

            thisvals = mean(suppression, 2, 'omitnan'); % average across stimuli

        case 'singletrial lifetime kurtosis'
            checkn_equals(n_equals, 'cells');
            thisvals = calculateLifetimeKurtosisSingleTrials(data);
        case 'stimuli lifetime kurtosis'
            checkn_equals(n_equals, 'cells');
            checkexists(stim_types, 'StimTypes');
            thisvals = calculateLifetimeKurtosisStimuli(data,stim_types);

        % metrics that return a matrix of size [time/units x trials]
        case 'max intensity'
            thisvals = squeeze(max(data,[],1,'omitmissing'));
        case 'avg intensity'
            thisvals = squeeze(mean(data,1,'omitmissing'));
        case 'variance'
            thisvals = squeeze(std(data,[],1,'omitmissing')).^2;
            
        % metrics that return a matrix of size [trials x trials]
        case 'common active units' % # TODO: common active units

        case 'mahalanobis distance' % # TODO: mahalanobis distance

        otherwise
            error('Unknown metric: %s', metric);
    end

    vals = thisvals; % Store the computed values

end

function suppression = getSuppressionScores(meanActivity)
    % Get the dimensions of the meanActivity matrix
    [numNeurons, numStims, numReps] = size(meanActivity);

    meanActivity = nanzscore(meanActivity, [], 3); % z-score across reps
    meanActivity = diff(meanActivity, 1, 3); % calculate difference across repetitions

    suppression = sum(meanActivity, 3, 'omitnan'); % sum over repetitions
    % [numNeurons x numStims]

end

function sparseness = calculateTuningSelectivity(tuningCurves)
    % Get the dimensions of the inputs
    [numNeurons, numStims, numRepetitions] = size(tuningCurves);

    % Initialize output vector: vertical vector with one [0->1] tuning
    % selectivity entry per cell
    sparseness = nan(numNeurons,numRepetitions);

    for i_rep = 1:numRepetitions
        meanActivity = tuningCurves(:,:,i_rep); % cells x stimuli

        % Calculate and store tuning selectivity for each neuron
        for i_cell = 1:numNeurons
            numerator = sum(meanActivity(i_cell,:)./numStims,'omitnan');
            denominator = sum((meanActivity(i_cell,:)).^2./numStims,'omitnan');
            sparseness(i_cell,i_rep) = 1 - (numerator^2 / denominator);
        end
    end
end

function sparseness = calculateLifetimeKurtosisSingleTrials(activityTraces)

     % Get the dimensions of the activityTraces matrix
    [timePoints, numNeurons, numTrials] = size(activityTraces);

    % Initialize output vector: vertical vector with one [0->1] sparseness
    % entry per trial
    sparseness = nan(numNeurons,1);

    % Calculate mean activity level for each neuron during each trial
    meanActivity = squeeze(nanmean(activityTraces, 1)); % cells x trials
    
    % Calculate and store lifetime kurtosis for each neuron
    sparseness = (nansum((zscore(meanActivity,[],2)).^4, 2) ./ numTrials ) -3;
end

function sparseness = calculateLifetimeKurtosisStimuli(activityTraces,stim_type)
    % Get the dimensions of the inputs
    [timePoints, numNeurons, numTrials] = size(activityTraces);
    stims = unique(stim_type);
    numStims = numel(stims);

    % Initialize output vector: vertical vector with one [0->1] tuning
    % selectivity entry per cell
    sparseness = nan(numNeurons,1);

    % Calculate mean activity level for each neuron for each stimulus
    allMeanActivity = squeeze(mean(activityTraces, 1,"omitmissing")); % cells x trials
    meanActivity = nan(numNeurons,numStims); % cells x stimuli
    for i_stim = 1:numStims
        thisstim = stims(i_stim);
        idx = ismember(stim_type,thisstim);
        
        tmp = mean(allMeanActivity(:,idx),2,'omitmissing');
        meanActivity(:,i_stim) = tmp;
    end
    
    % Calculate and store lifetime kurtosis for each neuron
    sparseness = (sum((zscore(meanActivity,[],2)).^4, 2,'omitmissing') ./ numStims ) -3;
end

function [tuningCurves, avgTuningCurves] = getTuningCurves(activityTraces, stim_type)
    % Returns tuning curves: [units x stimuli x repetitions]
    if nargin < 2
        error('getTuningCurves requires activityTraces and stim_type.');
    end

    [~, numUnits, numTrials] = size(activityTraces);
    stims = unique(stim_type);
    numStims = numel(stims);

    % Find repetitions per stimulus
    repCounts = arrayfun(@(s) sum(ismember(stim_type,s)), stims);
    maxReps = max(repCounts);

    tuningCurves = nan(numUnits, numStims, maxReps);

    for i_stim = 1:numStims
        stim = stims(i_stim);
        idx = find(ismember(stim_type,stim));
        nReps = numel(idx);
        for r = 1:nReps
            % Mean activity for each unit in this repetition
            tuningCurves(:, i_stim, r) = squeeze(mean(activityTraces(:, :, idx(r)), 1, 'omitmissing'));
        end
    end

    % Average tuning curves across repetitions
    avgTuningCurves = squeeze(mean(tuningCurves, 3, 'omitmissing'));
end

function checkn_equals(n_equals, expected)
    if ~strcmpi(n_equals, expected)
        error('Expected n_equals to be "%s", but got "%s".', expected, n_equals);
    end
end

function checkexists(var, name)
    if isempty(var)
        error('Expected %s to be provided.', name);
    end
end

% intermediate level plotter : dynamics
% 
% mode 'curve' : loop processor over a sliding window and plot values on a curve
% mode 'windows_size' : loop processor over a sliding window and over window sizes, plot imagesc onto provided axes
% mode 'repetitions' : loop processor over a sliding window and over repetitions, plot imagesc onto provided axes

% high level plotter : dynamics
% 
% compare across subject groups and stimulus groups, plot imagesc
% use colorbar, standard clim, log c?
