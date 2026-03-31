function vals = extractActivityMetric(events, metric, n_equals, varargin)
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
stim_types       = [];   % cell/array of stimulus labels, one per trial
baseline_mu      = [];   % [nUnits x 1] external baseline mean  (for responsiveness)
baseline_std     = [];   % [nUnits x 1] external baseline std   (for responsiveness)
pre_event_frames = [];   % integer: frames before event onset    (for motion tuning)
full_data        = [];   % [T x nUnits x n_trials] full trace    (for motion tuning event-free baseline)
event_onsets     = [];   % struct: .trial_idx, .time_s           (for motion tuning event-free baseline)
event_offsets    = [];   % struct: .trial_idx, .time_s           (for motion tuning event-free baseline)
fs               = [];   % scalar frame rate (Hz)                (for motion tuning event-free baseline)
ps_lim           = [];   % [t0, t1] peri-event window            (for motion tuning event-free baseline)
offset_padding_s = 4;    % seconds of padding after each offset  (for motion tuning event-free baseline)
if ~isempty(varargin)
    for i = 1:2:numel(varargin)
        switch lower(varargin{i})
            case 'stimtypes';       stim_types       = varargin{i+1};
            case 'baselinemu';      baseline_mu      = varargin{i+1};
            case 'baselinestd';     baseline_std     = varargin{i+1};
            case 'preeventframes';  pre_event_frames = varargin{i+1};
            case 'fulldata';        full_data        = varargin{i+1};
            case 'eventonsets';     event_onsets     = varargin{i+1};
            case 'eventoffsets';    event_offsets    = varargin{i+1};
            case 'framerate';       fs               = varargin{i+1};
            case 'pslim';           ps_lim           = varargin{i+1};
            case 'offsetpadding';   offset_padding_s = varargin{i+1};
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
        
        thisvals = calculateNormalizedSparseness(data, 1);
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
        thisvals = calculateNormalizedSparseness(meanActivity, 2); % lifetime sparseness using average activity across repetitions


    % metrics that return a column vector of size [units x 1]
    case 'selectivity of tuning'
        checkn_equals(n_equals, 'cells');
        checkexists(stim_types, 'StimTypes');

        % Calculate mean activity level for each neuron for each stimulus
        [~, meanActivity] = getTuningCurves(data, stim_types);

        % Calculate and store tuning selectivity for each neuron
        thisvals = calculateNormalizedSparseness(meanActivity, 2); % lifetime sparseness using average tuning curves

    case 'stability of tuning'
        checkn_equals(n_equals, 'cells');
        checkexists(stim_types, 'StimTypes');

        % Calculate mean activity level for each neuron for each stimulus
        meanActivity = getTuningCurves(data, stim_types); % [numNeurons x numStims x numReps]

        % Calculate average distance between tuning curves of all repetitions for each neuron
        avgDistance = nan(nUnits, 1);
        for i_unit = 1:nUnits
            tuningReps = squeeze(meanActivity(i_unit, :, :)); % [stimuli x repetitions]
            distances = pdist(tuningReps', 'correlation'); % pairwise distances between repetitions (sensitive to NaNs)
            
            avgDistance(i_unit) = mean(distances, 'omitnan'); % average distance
        end

        % Calculate and store tuning selectivity for each neuron
        thisvals = 1 - avgDistance; % [-1 -> 1], higher = more stable, 0 = uncorrelated, -1 = anticorrelated

    case 'stimulus shuffled selectivity of tuning'
        checkn_equals(n_equals, 'cells');
        checkexists(stim_types, 'StimTypes');

        % Calculate selectivity of stimulus-shuffled data
        nshuffles = 100;
        shuffledSelectivity = nan(nUnits, nshuffles);
        for i_shuffle = 1:nshuffles
            shuffledStimTypes = stim_types(randperm(length(stim_types)));
            [~, shuffledMeanActivity] = getTuningCurves(data, shuffledStimTypes);
            shuffledSelectivity(:, i_shuffle) = calculateNormalizedSparseness(shuffledMeanActivity, 2); % lifetime sparseness using average tuning curves
        end

        thisvals = shuffledSelectivity; % return

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
        % suppression = suppression - mean(suppression(2:nStims), 2, 'omitnan'); % subtract mean suppression score of the non-maximum stimuli
        
        thisvals = suppression(:,1) ./ std(suppression(:,2:nStims), 0, 2, 'omitnan'); % divide by STD of the non-maximum stimuli

    case 'general suppression score'
        % answers the question: how much is the activity of a neuron 
        % suppressed over repetitions on average?
        checkn_equals(n_equals, 'cells');
        checkexists(stim_types, 'StimTypes');
        
        % Calculate mean activity level for each neuron for each stimulus
        meanActivity = getTuningCurves(data, stim_types); % [numNeurons x numStims x numReps]

        suppression = getSuppressionScores(meanActivity); % [numNeurons x numStims]

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

    % -----------------------------------------------------------------------
    % responsiveness metrics — return [units x 1]
    %
    % Both metrics express unit responses as z-scores relative to a baseline:
    %   z = (mean_response - baseline_mu) / baseline_std
    %
    % 'odor responsiveness'
    %   Response window  : the events as provided (set interval to odor window
    %                      before calling).
    %   Baseline         : must be supplied by the caller via 'BaselineMu' and
    %                      'BaselineStd', computed from a full-session call
    %                      (interval = []).  See Figure1.m for reference.
    %
    % 'motion tuning'
    %   Response window  : frames after the event onset (t >= 0).
    %                      Requires 'PreEventFrames' to locate t = 0.
    %   Baseline (pref.) : subject-wide event-free intervals — all frames
    %                      across all trials that fall outside any motion
    %                      event window [onset, offset + OffsetPadding (4 s)].
    %                      Requires 'FullData', 'EventOnsets', 'EventOffsets',
    %                      'FrameRate'.  'PsLim' must match the window used
    %                      when detecting events.  'OffsetPadding' overrides
    %                      the default 4 s pad.
    %   Baseline (fallb.): pre-onset frames in the snippets if FullData is
    %                      absent but PreEventFrames is set.
    %   Baseline (ext.)  : 'BaselineMu'/'BaselineStd' if neither of the
    %                      above is available.
    % -----------------------------------------------------------------------

    case 'odor responsiveness'
        checkn_equals(n_equals, 'cells');
        thisvals = computeResponsiveness(data, baseline_mu, baseline_std);

    case 'motion tuning'
        checkn_equals(n_equals, 'cells');
        thisvals = computeMotionTuning(data, pre_event_frames, baseline_mu, baseline_std, ...
            full_data, event_onsets, event_offsets, fs, ps_lim, offset_padding_s);

    otherwise
        error('Unknown metric: %s', metric);
end

vals = thisvals; % Store the computed values

end

% =========================================================================
%  RESPONSIVENESS HELPERS
% =========================================================================

function vals = computeResponsiveness(data, baseline_mu, baseline_std)
% COMPUTERESPONSIVENESS  Per-unit responsiveness z-score.
%
%   Measures how many baseline standard deviations each unit's mean
%   activity during the supplied window exceeds the baseline mean.
%   Mirrors the calculation in Figure1.m:
%
%     resp = (mean_activity_in_window - mu_full) / std_full
%
%   INPUTS
%     data          [T x nUnits x nTrials]  peri-event activity
%     baseline_mu   [nUnits x 1]  session-wide mean, from a full-session
%                                 avg-intensity call (interval = [])
%     baseline_std  [nUnits x 1]  session-wide mean temporal std, from a
%                                 full-session variance call (interval = [])
%
%   OUTPUT
%     vals  [nUnits x 1]  responsiveness z-score per unit

    checkexists(baseline_mu,  'BaselineMu');
    checkexists(baseline_std, 'BaselineStd');

    % Mean over time → [nUnits x nTrials]; then mean over trials → [nUnits x 1]
    response = mean(squeeze(mean(data, 1, 'omitmissing')), 2, 'omitmissing');

    vals = (response(:) - baseline_mu(:)) ./ baseline_std(:);
end


function vals = computeMotionTuning(data, pre_event_frames, baseline_mu, baseline_std, ...
        full_data, event_onsets, event_offsets, fs, ps_lim, offset_padding_s)
% COMPUTEMOTIONTUNING  Per-unit motion-onset tuning z-score.
%
%   Measures each unit's response to tail-motion onset relative to a
%   subject-wide event-free baseline.
%
%   BASELINE PRIORITY
%     1. Event-free baseline (preferred): pool all frames across all trials
%        that fall outside any motion-event window [onset, offset + padding_s].
%        Requires 'FullData', 'EventOnsets', 'EventOffsets', 'FrameRate'.
%     2. Pre-event baseline (fallback): frames before t = 0 in the peri-event
%        snippets, as before.  Used when FullData is not supplied.
%     3. External baseline: 'BaselineMu' / 'BaselineStd' passed explicitly.
%
%   INPUTS
%     data              [T x nUnits x nEvents]  peri-onset snippets
%     pre_event_frames  integer  frames before onset in each snippet (t < 0).
%                       Still required to isolate the post-onset response
%                       regardless of which baseline path is taken.
%     baseline_mu       [nUnits x 1]  external fallback baseline mean
%     baseline_std      [nUnits x 1]  external fallback baseline std
%     full_data         [T_full x nUnits x n_trials]  full subject trace
%     event_onsets      struct  .trial_idx [n x 1], .time_s [n x 1]
%     event_offsets     struct  .trial_idx [n x 1], .time_s [n x 1]
%     fs                scalar  frame rate (Hz)
%     ps_lim            [t0, t1] or []  peri-event window boundaries
%     offset_padding_s  scalar  seconds to pad after each event offset
%
%   OUTPUT
%     vals  [nUnits x 1]  motion tuning z-score per unit

    % ------------------------------------------------------------------
    % 1. Determine baseline
    % ------------------------------------------------------------------
    use_event_free = ~isempty(full_data) && ~isempty(event_onsets) && ...
                     ~isempty(event_offsets) && ~isempty(fs);

    if use_event_free
        [T_full, nUnits, n_trials] = size(full_data);
        padding_fr = round(offset_padding_s * fs);

        % Build a logical free-frame mask: true where no event is active.
        free_mask = true(T_full, n_trials);

        for j = 1:n_trials
            on_idx = find(event_onsets.trial_idx == j);
            of_idx = find(event_offsets.trial_idx == j);

            if isempty(on_idx); continue; end

            % Convert times → frame indices (1-based)
            if isempty(ps_lim)
                on_fr = round(event_onsets.time_s(on_idx) * fs) + 1;
                of_fr = round(event_offsets.time_s(of_idx) * fs) + 1;
            else
                on_fr = round((event_onsets.time_s(on_idx)  - ps_lim(1)) * fs) + 1;
                of_fr = round((event_offsets.time_s(of_idx) - ps_lim(1)) * fs) + 1;
            end

            on_fr = sort(on_fr(:));
            of_fr = sort(of_fr(:));

            % Greedy pairing: each onset → first offset at or after it
            for k = 1:numel(on_fr)
                f_on   = on_fr(k);
                match  = of_fr(of_fr >= f_on);
                if isempty(match)
                    f_off = T_full;   % no offset found → extend to end of trial
                else
                    f_off = match(1);
                end
                f_start = max(1, f_on);
                f_end   = min(T_full, f_off + padding_fr);
                free_mask(f_start:f_end, j) = false;
            end
        end

        % Pool free frames: [nUnits x n_free_frames]
        free_flat = reshape(permute(full_data, [2, 1, 3]), nUnits, []);
        mask_flat = reshape(free_mask, 1, []);
        free_act  = free_flat(:, mask_flat);   % [nUnits x n_free_frames]

        bl_mu  = mean(free_act, 2, 'omitmissing');   % [nUnits x 1]
        bl_std = std(free_act,  0, 2, 'omitmissing'); % [nUnits x 1]

    elseif ~isempty(pre_event_frames) && pre_event_frames > 0
        % Fallback: pre-onset portion of the peri-event snippets
        pre   = data(1:pre_event_frames, :, :);   % [pre_fr x nUnits x nEvents]
        bl_mu  = mean(squeeze(mean(pre, 1, 'omitmissing')), 2, 'omitmissing');
        bl_std = mean(squeeze(std(pre,  0, 1, 'omitmissing')), 2, 'omitmissing');

    else
        % External baseline
        checkexists(baseline_mu,  'BaselineMu  (required when FullData and PreEventFrames are not provided)');
        checkexists(baseline_std, 'BaselineStd (required when FullData and PreEventFrames are not provided)');
        bl_mu  = baseline_mu(:);
        bl_std = baseline_std(:);
    end

    % ------------------------------------------------------------------
    % 2. Compute response (post-onset mean)
    % ------------------------------------------------------------------
    if ~isempty(pre_event_frames) && pre_event_frames > 0
        post = data(pre_event_frames+1:end, :, :);   % [post_fr x nUnits x nEvents]
    else
        post = data;
    end
    response = mean(squeeze(mean(post, 1, 'omitmissing')), 2, 'omitmissing'); % [nUnits x 1]

    % Guard against silent units (zero std → undefined tuning)
    bl_std(bl_std == 0) = NaN;

    vals = (response(:) - bl_mu(:)) ./ bl_std(:);
end


function sparseness = calculateNormalizedSparseness(meanActivity, dim)
    arguments
        meanActivity double % [neurons x stimuli] or [neurons x trials]
        dim (1,1) double {mustBeMember(dim, [1,2])} = 1 % dimension to operate along
    end

    % Handle 3D input by iterating over slices
    if ndims(meanActivity) == 3
        [dim1, dim2, dim3] = size(meanActivity);
        remainingdim = dim2;
        if dim==2; remainingdim = dim1; end
        sparseness = nan(remainingdim, dim3);
        for i = 1:dim3
            sparseness(:, i) = calculateNormalizedSparseness(meanActivity(:, :, i), dim);
        end
        return;
    end

    % If operating along dimension 2, transpose the matrix
    if dim == 2
        meanActivity = meanActivity'; % [stimuli x neurons] or [trials x neurons]
    end

    % Get the dimensions of the meanActivity matrix
    [numVars, numSamples] = size(meanActivity);

    % Calculate normalized sparseness for each column
    sparseness = nan(numSamples, 1);
    for i = 1:numSamples
        numerator = sum(meanActivity(:, i) ./ numVars, 'omitnan');
        denominator = sum((meanActivity(:, i)).^2 ./ numVars, 'omitnan');
        sparseness(i) = (1 - numerator^2 / denominator) / (1 - 1 / numVars);
    end

end

function suppression = getSuppressionScores(meanActivity)
    % Get the dimensions of the meanActivity matrix
    [numNeurons, numStims, numReps] = size(meanActivity);

    meanActivity = nanzscore(meanActivity, [], 3); % z-score across reps
    meanActivity = diff(meanActivity, 1, 3); % calculate difference across repetitions

    suppression = sum(meanActivity, 3, 'omitnan'); % sum over repetitions
    
    suppression = -suppression; % positive = more suppression
    % [numNeurons x numStims]

end

function sparseness = calculateTuningSelectivity(tuningCurves) % DEPRECATED
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
