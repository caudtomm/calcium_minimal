function [results] = extractSignalAndNoise(traces3d, stimlabs, varargin)
    % extractSignalAndNoise: Extract signal, noise, and confidence from trial data
    % 
    % Inputs:
    %   traces3d: A time x units x trials matrix containing trial data
    %   stimlabs: A 1D cell array of trial labels, defining stimulus label of that trial
    %   varargin: Additional parameters (unused in this minimal version)
    % 
    % Outputs:
    %   signalMAT: A matrix with computed signal for each trial
    %   confidenceMAT: A matrix with computed confidence for each trial
    %   noiseMAT: A matrix with computed noise for each trial
    %   signaltonoisecorrelations: Correlation between signal and noise
    %   original_trial_idx: Original indices of trials

    % For simplicity, in this minimal implementation:
    % Signal is the mean across units
    % Noise is the standard deviation across units
    % Confidence is arbitrarily set to 1 for this minimal example

    % Determine the size of the traces3d matrix
    [T, U, trialCount] = size(traces3d);

    stims = unique(stimlabs);
    nstims = numel(stims);

    % initialize vars
    Avg_traces = nan(T,U,nstims);
    VBL_traces = nan(T,U,nstims);
    Avg_vect = nan(U,nstims);
    Residuals_traces = nan(T,U,trialCount);
    Residuals_timeavg = nan(U,trialCount);
    original_trial_idx = false(trialCount,nstims);

    % Compute
    for i_stim = 1:nstims
        thisstim = stims{i_stim};
        idx = squeeze(ismember(stimlabs,thisstim));
        trialsdonesofar = sum(original_trial_idx,"all");
        
        % Get one matrix per stimulus type (odor), containing the
        % frame-by-frame trace average. Dim: [time x units]
        thisavg_f2f = nanmean(traces3d(:,:,idx),3);

        % Get one matrix per stimulus type (odor), containing the
        % frame-by-frame trace variability. Dim: [time x units]
        thisvariability_f2f = std(traces3d(:,:,idx),[],3,'omitnan').^2; % inter-trial variance

        % Get one vector per stimulus type (odor), containing the
        % overall time-averaged activity. Dim: [units]
        thisavg_all = squeeze(nanmean(thisavg_f2f,1));

        % Compute the residuals ("noise") around the average for each
        % trial. Dim: [time x units x trials]
        thisresiduals = traces3d(:,:,idx) - thisavg_f2f;

        % Compute the time-average of the frame-by-frame residuals. Dim:
        % [units x trials]
        thisresiduals_timeavg = squeeze(nanmean(thisresiduals,1));

        % Store in common matrices
        Avg_traces(:,:,i_stim) = thisavg_f2f;
        Avg_vect(:,i_stim) = thisavg_all;
        VBL_traces(:,:,i_stim) = thisvariability_f2f;
        Residuals_traces(:,:,trialsdonesofar+[1:sum(idx)]) = thisresiduals;
        Residuals_timeavg(:,trialsdonesofar+[1:sum(idx)]) = thisresiduals_timeavg;
        original_trial_idx(:,i_stim) = idx;
    end

    % get index list of trials, arranged by stimulus type and in
    % chronological order within each stimulus type
    idx_by_stim_type = [];
    for i = 1:nstims
        idx_by_stim_type = [idx_by_stim_type ; ...
            squeeze(find(original_trial_idx(:,i)))];
    end

    % Store for return
    results.signalMAT = Avg_traces;
    results.avg_signal_vectors = Avg_vect;
    results.confidenceMAT = 1./VBL_traces;
    results.noiseMAT = Residuals_traces;
    results.noiseMAT_timeavg = Residuals_timeavg;
    results.original_trial_idx = original_trial_idx;
    results.idx_by_stim_type = idx_by_stim_type;
    results.stim_order = stims;


    %% Plotting
    
    



end
