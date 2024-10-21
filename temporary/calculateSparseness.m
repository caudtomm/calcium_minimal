function sparseness = calculateSparseness(activityTraces,method,varargin)
    % Input:
    % - activityTraces: Time x Neurons x Trials matrix of neuronal activity traces

    %% parse inputs

    % Set default values for the optional arguments

    [timePoints, numNeurons, numTrials] = size(activityTraces);
    stim_type = repelem({'default'},numTrials);
    
    % Parse the varargin input and update the corresponding arguments
    numArgs = numel(varargin);
    if mod(numArgs, 2) ~= 0
        error('Invalid number of arguments. Please provide name=value pairs.');
    end
    
    for i = 1:2:numArgs
        name = varargin{i};
        value = varargin{i+1};
        
        switch name
            case 'StimTypes'
                stim_type = value;
            otherwise
                error('Invalid argument name: %s', name);
        end
    end

    %% initialize output
    sparseness = [];

    switch method
        case 'population'
            sparseness = calculatePopulationSparseness(activityTraces);
        case 'selectivity'
            sparseness = calculateTuningSelectivity(activityTraces,stim_type);
        case 'kurtosis-singletrials'
            sparseness = calculateLifetimeKurtosisSingleTrials(activityTraces);
        case 'kurtosis-stimuli'
            sparseness = calculateLifetimeKurtosisStimuli(activityTraces,stim_type);
        otherwise
            warning('sparseness method unknown! returned blank')
            return
    end

   
end


function sparseness = calculatePopulationSparseness(activityTraces)
    % Get the dimensions of the activityTraces matrix
    [timePoints, numNeurons, numTrials] = size(activityTraces);

    % Initialize output vector: vertical vector with one [0->1] sparseness
    % entry per trial
    sparseness = nan(numTrials,1);

    % Calculate mean activity level for each neuron during each trial
    meanActivity = squeeze(nanmean(activityTraces, 1)); % cells x trials
    
    % Calculate and store population sparseness for each trial
    for i_trial = 1:numTrials
        numerator = sum(meanActivity(:,i_trial)./numNeurons,'omitnan');
        denominator = sum((meanActivity(:,i_trial)).^2./numNeurons,'omitnan');
        sparseness(i_trial) = 1 - (numerator^2 / denominator);
    end
end

function sparseness = calculateTuningSelectivity(activityTraces,stim_type)
    % Get the dimensions of the inputs
    [timePoints, numNeurons, numTrials] = size(activityTraces);
    stims = unique(stim_type);
    numStims = numel(stims);

    % Initialize output vector: vertical vector with one [0->1] tuning
    % selectivity entry per cell
    sparseness = nan(numNeurons,1);

    % Calculate mean activity level for each neuron for each stimulus
    allMeanActivity = squeeze(nanmean(activityTraces, 1)); % cells x trials
    meanActivity = nan(numNeurons,numStims); % cells x stimuli
    for i_stim = 1:numStims
        thisstim = stims(i_stim);
        idx = ismember(stim_type,thisstim);
        
        tmp = nanmean(allMeanActivity(:,idx),2);
        meanActivity(:,i_stim) = tmp;
    end
    
    % Calculate and store tuning selectivity for each neuron
    for i_cell = 1:numNeurons
        numerator = sum(meanActivity(i_cell,:)./numStims,'omitnan');
        denominator = sum((meanActivity(i_cell,:)).^2./numStims,'omitnan');
        sparseness(i_cell) = 1 - (numerator^2 / denominator);
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
    allMeanActivity = squeeze(nanmean(activityTraces, 1)); % cells x trials
    meanActivity = nan(numNeurons,numStims); % cells x stimuli
    for i_stim = 1:numStims
        thisstim = stims(i_stim);
        idx = ismember(stim_type,thisstim);
        
        tmp = nanmean(allMeanActivity(:,idx),2);
        meanActivity(:,i_stim) = tmp;
    end
    
    % Calculate and store lifetime kurtosis for each neuron
    sparseness = (nansum((zscore(meanActivity,[],2)).^4, 2) ./ numStims ) -3;
end





