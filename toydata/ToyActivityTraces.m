classdef ToyActivityTraces
% ToyActivityTraces  ActivityTraces-compatible adapter for toy model data.
%
%   Wraps a single ToyDataGenerator .mat file into an object that exposes
%   the same public interface as ActivityTraces, so it can be stored in
%   Experiment.traces{} and consumed by ExperimentViewer, DataFilter,
%   ModeSelector, and GCMC_Analysis without modification.
%
%   Dimension conventions follow ActivityTraces:
%     pSpike / dFoverF : [T x N x T_trials]  where T_trials = numel(stim_idx)
%     stim_series      : table with columns trialnum, odor_channel,
%                        frame_onset, frame_offset, stimulus, rep
%
%   Trial ordering follows params.stim_idx exactly.
%   Default stim_idx (empty): repelem(1:K, R)  (odor-major, R reps each).
%
%   Usage:
%     s   = load('toydata_hyp_...mat');
%     at  = ToyActivityTraces(s, 'trained');
%     at.pSpike          % [T x N x T_trials]
%     at.stim_series     % table
%     at.dFoverF_good    % [T x N x T_trials]  (dependent, mirrors pSpike)
%
%   See also: ToyDataLoader, ToyDataGenerator, ActivityTraces

    % ----------------------------------------------------------------- %
    % Core interface properties (mirror ActivityTraces naming exactly)
    % ----------------------------------------------------------------- %
    properties

        % Primary data: inferred firing rates
        pSpike          double  % [T x N x ntrials]  firing rates (Hz)
        dFoverF         double  % [T x N x ntrials]  = pSpike (for dFoverF pipeline)

        % Experimental metadata
        stim_series     table   % trialnum, odor_channel, frame_onset, frame_offset, stimulus
        subject_group   char    % group label (e.g. 'trained', 'naive', 'toy')
        framerate       double  % sampling rate (Hz)
        t               double  % [1 x T]  time axis (seconds)

        % Dimensions
        ntrials         double  % K * R
        N               double  % number of neurons
        L               double  % frames per trial
        T               double  % seconds per trial

        % Quality control (empty for toy data — all neurons / trials are good)
        goodNeuron_IDs  double  % [1 x N]  1:N
        ROImap          double  % [1 x N]  synthetic identity map
        badtrials       double  % []
        badperiods              % []

        % Locations (empty stub — toy data has no file-backed Subject)
        subject_locations Locations

    end

    % ----------------------------------------------------------------- %
    % Dependent properties (match ActivityTraces interface exactly)
    % ----------------------------------------------------------------- %
    properties (Dependent)
        dFoverF_good    double  % [T x N x ntrials]  = dFoverF(:, goodNeuron_IDs, :)
    end

    % ----------------------------------------------------------------- %
    % Toy-specific extras (ignored by the pipeline, useful for inspection)
    % ----------------------------------------------------------------- %
    properties
        toy_params      struct  % full params struct from the .mat file
        toy_metadata    struct  % metadata struct from the .mat file
        toy_geometry    struct  % geometry struct (e, n, d, U)
        toy_patterns    struct  % patterns struct (mu, novelty_weights, etc.)
    end

    % --------------------------------------------------------------------- %
    methods

        function obj = ToyActivityTraces(s, group)
        % TOYACTIVITYTRACES  Construct adapter from a loaded .mat struct.
        %
        %   at = ToyActivityTraces(s, group)
        %
        %   s     : struct loaded from ToyDataGenerator .mat file
        %   group : char  group label to assign (e.g. 'trained', 'toy_hyp')

            if nargin < 2 || isempty(group)
                group = 'toy';
            end

            meta = s.metadata;                % struct from build_metadata
            fr   = s.firing_rates;            % [N x T x T_trials]

            N_neurons = meta.N;
            T_frames  = meta.T;
            T_trials  = meta.T_trials;

            % Permute: [N x T x T_trials] → [T x N x T_trials]
            fr_perm = permute(fr, [2, 1, 3]);

            % Core data
            obj.pSpike  = fr_perm;
            obj.dFoverF = fr_perm;

            % Temporal properties
            obj.framerate = meta.fs;
            obj.t         = meta.time;        % [1 x T]
            obj.L         = T_frames;
            obj.T         = T_frames / meta.fs;

            % Dimensions
            obj.ntrials = T_trials;
            obj.N       = N_neurons;

            % Quality control (all neurons good, no bad trials)
            obj.goodNeuron_IDs = 1 : N_neurons;
            obj.ROImap         = reshape(1 : N_neurons, 1, N_neurons);
            obj.badtrials      = [];
            obj.badperiods     = [];

            % Group and locations
            obj.subject_group     = group;
            obj.subject_locations = Locations();

            % Stimulus table
            obj.stim_series = build_stim_series(meta);

            % Extras for inspection
            obj.toy_params   = s.params;
            obj.toy_metadata = meta;
            if isfield(s, 'geometry'); obj.toy_geometry = s.geometry; end
            if isfield(s, 'patterns'); obj.toy_patterns = s.patterns; end
        end

        % -- Dependent getter ----------------------------------------- %

        function val = get.dFoverF_good(obj)
        % Mirrors ActivityTraces: dFoverF(:, goodNeuron_IDs, :)
            val = obj.dFoverF(:, obj.goodNeuron_IDs, :);
        end

    end

end

% ========================================================================= %
% Local helper
% ========================================================================= %

function stim_series = build_stim_series(meta)
% BUILD_STIM_SERIES  Construct a stim_series table matching ActivityTraces format.
%
%   Uses meta.stim_idx for the exact trial presentation order.
%   Columns: trialnum, odor_channel, frame_onset, frame_offset, stimulus, rep

    idx          = meta.stim_idx(:);   % [T_trials x 1]
    T_trials     = numel(idx);
    trialnum     = (1 : T_trials)';
    odor_channel = idx;
    frame_onset  = repmat(meta.odor_frames(1), T_trials, 1);
    frame_offset = repmat(meta.odor_frames(2), T_trials, 1);
    stimulus     = cell(T_trials, 1);
    rep          = zeros(T_trials, 1);

    rep_count = zeros(1, meta.K);
    for t = 1:T_trials
        k              = idx(t);
        rep_count(k)   = rep_count(k) + 1;
        rep(t)         = rep_count(k);
        stimulus{t}    = meta.stimulus_names{k};
    end

    stim_series = table(trialnum, odor_channel, frame_onset, frame_offset, stimulus, rep, ...
                        'VariableNames', ...
                        {'trialnum', 'odor_channel', 'frame_onset', 'frame_offset', 'stimulus', 'rep'});
end
