classdef ToyParams
% ToyParams  Parameter container for the funnel generative model.
%
%   Encodes all parameters needed to generate synthetic calcium imaging
%   data with a specified geometric structure. Odor representations live on
%   a "funnel" global manifold: novelty-driven decay toward the origin
%   along odor-specific axes, plus an optional shared drift component.
%
%   Null model:  B = 0, gamma = 0, lambda_A = 0, C = 0
%                (or use ToyParams.makeNull())
%
%   Usage:
%     p = ToyParams()                              % all defaults
%     p = ToyParams('tau', 2, 'rho', 0.3)         % name-value pairs
%     p = ToyParams(struct('tau', 2, 'rho', 0.3)) % from struct
%     p.validate()                                 % check consistency
%
%   See also: ToyDataGenerator, ToyDataLoader

    properties

        % ----------------------------------------------------------------
        % Population dimensions
        % ----------------------------------------------------------------

        N double = 200   % number of neurons
        K double = 6     % number of odors
        R double = 5     % repetitions per odor used to build the default stim_idx

        % stim_idx: [1 x T_trials] presentation order of odors (values in 1..K).
        %           Defines the full trial sequence; e.g. [1 2 3 1 2 3] means
        %           odor A, B, C, A, B, C in order.
        %           [] = default: repelem(1:K, R)  (odor-major, R reps each).
        %           Different odors may appear different numbers of times.
        stim_idx double = []

        % ----------------------------------------------------------------
        % Manifold geometry
        % ----------------------------------------------------------------

        % D: dimensionality of the identity subspace span(U).
        %    Determines how many directions are available for identity axes,
        %    and how many novelty / drift axes can be orthogonal to span(U).
        %    Pairwise correlations of e_k are controlled independently by rho_e.
        %    Constraints: K <= D <= N.
        D double = 30

        % A: identity amplitude (Hz), uniform across odors.
        %    Sets the magnitude of the odor-specific component at rep 1.
        A double = 5

        % lambda_A: fractional growth of the identity amplitude across reps.
        %           identity component = A * (1 + lambda_A*(1-exp(-(r-1)/tau))) * e_k
        %           0 = flat (default, recovers original model);
        %           1 = amplitude doubles asymptotically.
        %           Uses time constant tau_A.
        lambda_A double = 0

        % B: initial novelty amplitude (Hz) at repetition 1.
        %    Decays exponentially across repetitions with time constant tau.
        %    Set to 0 for the null model (no novelty-driven attenuation).
        B double = 8

        % tau: novelty decay time constant (units: repetitions).
        %      novelty component = B * exp(-(r-1) / tau).
        tau double = 1.5

        % tau_A: identity growth time constant (units: repetitions).
        %        Governs both lambda_A and eta saturation:
        %          identity scale = 1 + lambda_A*(1 - exp(-(r-1)/tau_A))
        %          rho_e_eff(r)   = rho_e + eta*(1 - exp(-(r-1)/tau_A))
        tau_A double = 1.5

        % gamma: directed drift strength (Hz per repetition) along the
        %        odor-specific drift axis d_k.
        %        Set to 0 for the null model (no drift).
        gamma double = 1

        % C: modulus (Hz) of the trial-specific baseline vector s_t.
        %    ||s_t|| = C exactly at every trial.
        %    0 = no baseline term (default).
        C double = 0

        % theta: baseline drift rate across trials (rotation speed).
        %        s_t = C*[cos(theta*pi/2*(t-1))*u_s + sin(theta*pi/2*(t-1))*v_s]
        %        0 = constant across trials (global offset).
        %        1 = consecutive trials are orthogonal (no correlation).
        %        Values in (0,1): linear angular drift, partial correlation.
        %        Values > 1: consecutive trials anti-correlated.
        theta double = 0

        % ----------------------------------------------------------------
        % Axis correlations
        % ----------------------------------------------------------------

        % alpha: inner product <e_k, n_k> within each odor.
        %        Controls alignment between identity and novelty axes.
        %        0 = orthogonal (default); range [-1, 1].
        alpha double = 0

        % rho: inner product <n_k, n_j> between novelty axes of different
        %      odors (k ~= j). 0 = uncorrelated (default); range [-1, 1].
        %      Must be >= -1/(K-1) to keep the Gram matrix PSD.
        rho double = 0

        % rho_e: starting inner product <e_k, e_j> between identity axes of
        %        different odors (k ~= j) at rep 1.
        %        0 = orthogonal (default); range [-1, 1].
        %        Must be >= -1/(K-1) to keep the Gram matrix PSD.
        %        Independent of D (which sets subspace dimensionality only).
        rho_e double = 0

        % eta: gain of the exponential saturation of rho_e across repetitions.
        %      rho_e_eff(r) = rho_e + eta * (1 - exp(-(r-1)/tau_A))
        %      0    = constant rho_e (default).
        %      eta > 0: identity axes grow more correlated over reps.
        %      eta < 0: identity axes grow less correlated (or more
        %               anti-correlated) over reps; useful when rho_e > 0
        %               or to push decorrelated axes toward anti-correlation.
        %      Constraint: rho_e + eta in [-1/(K-1), 1].
        eta double = 0

        % rho_d: inner product <d_k, d_j> between drift directions of
        %        different odors (k ~= j). 0 = uncorrelated (default);
        %        range [-1/(K-1), 1].  rho_d = 1 recovers a single shared
        %        drift direction (original behaviour).
        rho_d double = 0

        % ----------------------------------------------------------------
        % Mixed selectivity
        % ----------------------------------------------------------------

        % rotation_mix: controls how the abstract geometry is aligned to
        %               the neuron basis.
        %               0 = identity-aligned (sparse, single-odor cells);
        %               1 = fully random rotation (dense mixed selectivity).
        %               Range [0, 1].
        rotation_mix double = 0

        % ----------------------------------------------------------------
        % Noise model
        % ----------------------------------------------------------------

        % fano_factor: trial-to-trial Fano factor (variance / mean of
        %              per-trial mean firing rate during the odor window).
        fano_factor double = 1.2

        % fr_mean: target mean firing rate across neurons (Hz).
        %          Sets the overall scale of baseline and odor-window activity.
        fr_mean double = 2

        % fr_cv: coefficient of variation of the across-neuron firing rate
        %        distribution (lognormal model). Controls neuron-to-neuron
        %        heterogeneity. fr_cv = std(FR) / mean(FR).
        fr_cv double = 1

        % noise_corr: pairwise noise correlation between neurons.
        %             0 = independent; values in [0, 1).
        noise_corr double = 0

        % noise_tau_frames: time constant of the exponential smoothing
        %                   kernel applied to the noise (frames).
        %                   Controls temporal autocorrelation within a trial.
        noise_tau_frames double = 8

        % ----------------------------------------------------------------
        % Temporal structure
        % ----------------------------------------------------------------

        fs           double = 7.67  % sampling rate (Hz)
        T            double = 1300  % total frames per trial
        t_odor_start double = 30    % odor window start (s)
        t_odor_end   double = 50    % odor window end (s)

        % ----------------------------------------------------------------
        % Metadata
        % ----------------------------------------------------------------

        % stimulus_names: {1 x K} cell of char labels.
        %                 Defaults to {'odor1', ..., 'odorK'} if empty.
        stimulus_names cell = {}

        % rng_seed: integer seed for reproducibility across replicas.
        rng_seed double = 42

    end

    properties (Dependent)
        % odor_frames: [start_frame, end_frame] of the odor window (1-indexed).
        odor_frames double
    end

    % --------------------------------------------------------------------- %
    methods

        function obj = ToyParams(varargin)
            if nargin == 0
                return
            end
            if nargin == 1 && isstruct(varargin{1})
                s = varargin{1};
                fn = fieldnames(s);
                for i = 1:numel(fn)
                    if isprop(obj, fn{i})
                        obj.(fn{i}) = s.(fn{i});
                    else
                        warning('ToyParams: ignoring unknown field "%s".', fn{i});
                    end
                end
            elseif mod(nargin, 2) == 0
                for i = 1:2:nargin
                    name = varargin{i};
                    val  = varargin{i+1};
                    if isprop(obj, char(name))
                        obj.(char(name)) = val;
                    else
                        error('ToyParams: unknown property "%s".', char(name));
                    end
                end
            else
                error('ToyParams: use ToyParams(), ToyParams(struct), or ToyParams(name, value, ...).');
            end
        end

        % -- Dependent getter ------------------------------------------- %

        function frames = get.odor_frames(obj)
            f_start = round(obj.t_odor_start * obj.fs) + 1;
            f_end   = min(round(obj.t_odor_end * obj.fs), obj.T);
            frames  = [f_start, f_end];
        end

        % -- Validation ------------------------------------------------- %

        function validate(obj)
        % VALIDATE  Check cross-parameter consistency. Errors on failure.
            assert(obj.D >= obj.K, ...
                'ToyParams: D (%d) must be >= K (%d).', obj.D, obj.K);
            assert(obj.D <= obj.N, ...
                'ToyParams: D (%d) cannot exceed N (%d).', obj.D, obj.N);

            assert(obj.alpha >= -1 && obj.alpha <= 1, ...
                'ToyParams: alpha must be in [-1, 1] (got %.3f).', obj.alpha);
            assert(obj.rho >= -1 && obj.rho <= 1, ...
                'ToyParams: rho must be in [-1, 1] (got %.3f).', obj.rho);
            assert(obj.rho_e >= -1 && obj.rho_e <= 1, ...
                'ToyParams: rho_e must be in [-1, 1] (got %.3f).', obj.rho_e);

            % Minimum rho / rho_e for PSD Gram matrix: 1 + (K-1)*rho >= 0
            rho_min = -1 / (obj.K - 1);
            assert(obj.rho >= rho_min, ...
                ['ToyParams: rho (%.3f) too negative for K=%d. ' ...
                 'Minimum valid rho is %.3f.'], obj.rho, obj.K, rho_min);
            assert(obj.rho_e >= rho_min, ...
                ['ToyParams: rho_e (%.3f) too negative for K=%d. ' ...
                 'Minimum valid rho_e is %.3f.'], obj.rho_e, obj.K, rho_min);

            rho_e_asym = obj.rho_e + obj.eta;   % asymptotic rho_e_eff as r -> inf
            assert(rho_e_asym <= 1, ...
                ['ToyParams: rho_e + eta (%.3f) exceeds 1. ' ...
                 'Asymptotic rho_e_eff must be <= 1.'], rho_e_asym);
            assert(rho_e_asym >= rho_min, ...
                ['ToyParams: rho_e + eta (%.3f) violates PSD bound for K=%d. ' ...
                 'Minimum asymptotic rho_e_eff is %.3f.'], rho_e_asym, obj.K, rho_min);

            assert(obj.rho_d >= -1 && obj.rho_d <= 1, ...
                'ToyParams: rho_d must be in [-1, 1] (got %.3f).', obj.rho_d);
            assert(obj.rho_d >= rho_min, ...
                ['ToyParams: rho_d (%.3f) too negative for K=%d. ' ...
                 'Minimum valid rho_d is %.3f.'], obj.rho_d, obj.K, rho_min);

            assert(obj.rotation_mix >= 0 && obj.rotation_mix <= 1, ...
                'ToyParams: rotation_mix must be in [0, 1] (got %.3f).', obj.rotation_mix);

            assert(obj.noise_corr >= 0 && obj.noise_corr < 1, ...
                'ToyParams: noise_corr must be in [0, 1) (got %.3f).', obj.noise_corr);

            assert(obj.fano_factor > 0, ...
                'ToyParams: fano_factor must be positive.');
            assert(obj.fr_mean > 0, ...
                'ToyParams: fr_mean must be positive.');
            assert(obj.fr_cv > 0, ...
                'ToyParams: fr_cv must be positive.');
            assert(obj.tau > 0, ...
                'ToyParams: tau must be positive.');
            assert(obj.tau_A > 0, ...
                'ToyParams: tau_A must be positive.');
            assert(obj.B >= 0, ...
                'ToyParams: B must be non-negative.');
            assert(obj.gamma >= 0, ...
                'ToyParams: gamma must be non-negative.');
            assert(obj.lambda_A >= 0, ...
                'ToyParams: lambda_A must be non-negative.');
            assert(obj.C >= 0, ...
                'ToyParams: C must be non-negative.');
            assert(obj.theta >= 0, ...
                'ToyParams: theta must be non-negative.');

            if ~isempty(obj.stim_idx)
                assert(all(obj.stim_idx == floor(obj.stim_idx)) && ...
                       all(obj.stim_idx >= 1) && all(obj.stim_idx <= obj.K), ...
                    'ToyParams: stim_idx values must be integers in 1..K (%d).', obj.K);
            end

            assert(obj.t_odor_end > obj.t_odor_start, ...
                'ToyParams: t_odor_end must be > t_odor_start.');
            assert(obj.odor_frames(2) <= obj.T, ...
                'ToyParams: odor window end (frame %d) exceeds T=%d.', ...
                obj.odor_frames(2), obj.T);

            if ~isempty(obj.stimulus_names)
                assert(numel(obj.stimulus_names) == obj.K, ...
                    'ToyParams: numel(stimulus_names) must equal K (%d).', obj.K);
            end
        end

        % -- Convenience ----------------------------------------------- %

        function idx = get_stim_idx(obj)
        % GET_STIM_IDX  Returns stim_idx, with default if empty.
        %   Default: odor-major ordering, R repetitions each: [1..1 2..2 ... K..K]
            if isempty(obj.stim_idx)
                idx = repelem(1:obj.K, obj.R);
            else
                idx = obj.stim_idx;
            end
        end

        function names = get_stimulus_names(obj)
        % GET_STIMULUS_NAMES  Returns stimulus_names, with defaults if empty.
            if isempty(obj.stimulus_names)
                names = arrayfun(@(k) sprintf('odor%d', k), 1:obj.K, ...
                                 'UniformOutput', false);
            else
                names = obj.stimulus_names;
            end
        end

        function s = toStruct(obj)
        % TOSTRUCT  Convert to plain struct for .mat serialization.
            mc   = metaclass(obj);
            props = mc.PropertyList;
            s = struct();
            for i = 1:numel(props)
                if ~props(i).Dependent
                    s.(props(i).Name) = obj.(props(i).Name);
                end
            end
        end

        function p = makeNull(obj)
        % MAKENULL  Returns a null-model ToyParams from an existing ToyParams object.
        %   Sets B=0, gamma=0, lambda_A=0, C=0 (all rep- and trial-dependent
        %   terms removed); stim_idx, noise, and geometry params are preserved.
            p          = obj;
            p.B        = 0;
            p.gamma    = 0;
            p.lambda_A = 0;
            p.C        = 0;
        end

    end

    % --------------------------------------------------------------------- %
    methods (Static)

        function p = fromStruct(s)
        % FROMSTRUCT  Reconstruct a ToyParams object from a saved struct.
            p = ToyParams(s);
        end

    end

end
