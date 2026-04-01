classdef NoiseModel
% NoiseModel  Generates firing rate traces from central patterns and noise.
%
%   Adds temporally correlated, spatially correlated noise to the central
%   patterns produced by CentralPatterns, yielding a full [N x T x K x R]
%   array of synthetic firing rates.
%
%   Noise model:
%     noise_i(t) = sigma_i * [ sqrt(c)*eta(t) + sqrt(1-c)*xi_i(t) ]
%
%   where eta(t) is a shared colored noise process, xi_i(t) are independent
%   colored noise processes, and c = params.noise_corr.  Both eta and xi are
%   causal exponential (IIR) filtered white noise.  This gives:
%     Corr(noise_i(t), noise_j(t)) = noise_corr   for i~=j
%     Var(noise_i(t)) = sigma_i^2
%
%   Noise amplitude (sigma_i) is calibrated so that the trial-to-trial
%   variance of the time-averaged response matches the target Fano factor
%   at baseline activity levels:
%     sigma_i = sqrt(fano_factor * baseline_i / ac_factor)
%
%   where ac_factor = Var(mean_{odor window}(noise)) / Var(noise(t)).
%   This correction accounts for temporal autocorrelation reducing the
%   effective variance of the time-averaged noise.
%
%   Note: effective Fano factor during odor responses is
%     Fano_eff(i,k,r) = fano_factor * baseline_i / (baseline_i + mu(i,k,r))
%   which is less than the target for activated neurons (biologically realistic).
%
%   Usage:
%     params = ToyParams();
%     geom   = GeometryBuilder(params);
%     cp     = CentralPatterns(params, geom);
%     nm     = NoiseModel(params, cp);
%     nm.firing_rates   % [N x T x K x R]
%     nm.report()
%
%   See also: ToyParams, GeometryBuilder, CentralPatterns, ToyDataGenerator

    properties (SetAccess = private)

        firing_rates    double  % [N x T x T_trials]  firing rates (Hz), clipped to >= 0
        baseline        double  % [N x 1]           per-neuron mean baseline rate (Hz)
        sigma_noise     double  % [N x 1]           per-neuron noise std (Hz)
        ac_factor       double  % scalar  Var(mean_T_odor(noise)) / Var(noise(t))

    end

    % --------------------------------------------------------------------- %
    methods

        function obj = NoiseModel(params, cp)
        % NOISEMODEL  Generate firing rate traces.
        %
        %   nm = NoiseModel(params, cp)
        %
        %   params : ToyParams
        %   cp     : CentralPatterns
        %
        %   Continues from the RNG state left by GeometryBuilder.

            % 1. Per-neuron baseline firing rates (lognormal)
            obj.baseline = sample_baselines(params);

            % 2. Autocorrelation correction factor for the odor window
            obj.ac_factor = compute_ac_factor(params);

            % 3. Per-neuron noise amplitude (Fano-matched to baseline)
            obj.sigma_noise = sqrt(params.fano_factor .* obj.baseline ./ obj.ac_factor);

            % 4. Generate all noise: [N x T x T_trials]
            T_trials = size(cp.mu, 2);
            noise = generate_noise(params, obj.sigma_noise, T_trials);

            % 5. Assemble signal + noise, clip to >= 0
            obj.firing_rates = assemble_traces(params, cp, obj.baseline, noise);
        end

        % -- Diagnostics ------------------------------------------------ %

        function report(obj, params)
        % REPORT  Print summary statistics of the generated firing rates.
        %
        %   nm.report(params)

            fr = obj.firing_rates;   % [N x T x T_trials]
            [N, ~, T_trials] = size(fr);

            f1 = params.odor_frames(1);
            f2 = params.odor_frames(2);

            % Mean over odor window per (neuron, trial): [N x T_trials]
            odor_mean = squeeze(mean(fr(:, f1:f2, :), 2));

            % Effective Fano: trial-to-trial var / mean across all trials
            trial_var  = var(odor_mean, 0, 2);   % [N x 1]
            trial_mean = mean(odor_mean, 2);      % [N x 1]
            active     = trial_mean > 0.1;
            fano_est   = trial_var(active) ./ trial_mean(active);

            fprintf('NoiseModel diagnostics\n');
            fprintf('  T_trials = %d\n', T_trials);
            fprintf('  Baseline firing rates (Hz):  mean = %.3f,  std = %.3f\n', ...
                    mean(obj.baseline), std(obj.baseline));
            fprintf('  Noise sigma (Hz):            mean = %.3f,  std = %.3f\n', ...
                    mean(obj.sigma_noise), std(obj.sigma_noise));
            fprintf('  Autocorrelation factor:      %.6f\n', obj.ac_factor);
            fprintf('  Effective Fano (active neurons, across trials):\n');
            fprintf('    mean = %.3f,  median = %.3f,  target = %.3f\n', ...
                    mean(fano_est), median(fano_est), params.fano_factor);
            fprintf('  Firing rate range: [%.3f, %.3f] Hz\n', ...
                    min(fr(:)), max(fr(:)));
        end

    end

end

% ========================================================================= %
% Local helper functions
% ========================================================================= %

function baseline = sample_baselines(params)
% SAMPLE_BASELINES  Draw per-neuron baseline firing rates from lognormal.
%
%   Lognormal parameters derived from target mean and CV:
%     sigma_ln^2 = log(1 + CV^2)
%     mu_ln      = log(fr_mean) - sigma_ln^2 / 2

    sigma_ln = sqrt(log(1 + params.fr_cv^2));
    mu_ln    = log(params.fr_mean) - sigma_ln^2 / 2;
    baseline = lognrnd(mu_ln, sigma_ln, params.N, 1);
end

% ------------------------------------------------------------------------- %

function ac = compute_ac_factor(params)
% COMPUTE_AC_FACTOR  Exact Var(mean_T(eps)) / Var(eps(t)) for exp-colored noise.
%
%   For a stationary process with autocorrelation R(lag) = sigma^2 * alpha^|lag|,
%   the variance of its time-average over T frames is:
%
%     Var(mean_T) = sigma^2 / T^2 * sum_{lag=0}^{T-1} 2*(T-lag)*alpha^lag - sigma^2/T
%                = sigma^2 * sum_{lag=0}^{T-1} (T-lag) * alpha^lag * 2 / T^2  - ...
%
%   Computed exactly as:
%     ac = (1/T^2) * sum_{lag=0}^{T-1} count(lag) * alpha^lag
%   where count(lag) = T for lag=0, 2*(T-lag) for lag>0.

    T_odor = diff(params.odor_frames) + 1;
    alpha  = exp(-1 / params.noise_tau_frames);
    lags   = 0 : (T_odor - 1);

    % Number of (t1, t2) pairs with |t1-t2| = lag
    counts        = T_odor - lags;          % lag=0: T pairs, lag=k: T-k pairs
    counts(2:end) = 2 * counts(2:end);      % symmetric: both +lag and -lag

    ac = sum(counts .* alpha .^ lags) / T_odor^2;
end

% ------------------------------------------------------------------------- %

function noise = generate_noise(params, sigma_noise, T_trials)
% GENERATE_NOISE  Returns [N x T x T_trials] spatially and temporally correlated noise.
%
%   Steps:
%     1. Draw independent white noise [N x T x T_trials] and shared [1 x T x T_trials]
%     2. Mix with factor model for spatial correlation: noise_corr
%     3. Apply causal IIR filter along time (dim 2) for temporal correlation
%     4. Normalize to unit variance (correct for filter's variance scaling)
%     5. Scale each neuron by sigma_noise(i)

    N = params.N;
    T = params.T;
    c = params.noise_corr;

    alpha   = exp(-1 / params.noise_tau_frames);
    iir_b   = 1 - alpha;
    iir_a   = [1, -alpha];
    iir_std = sqrt((1-alpha) / (1+alpha));

    x_indep  = randn(N, T, T_trials);
    x_shared = randn(1, T, T_trials);   % broadcast over neurons

    x_mixed   = sqrt(c) * x_shared + sqrt(1 - c) * x_indep;   % [N x T x T_trials]
    x_colored = filter(iir_b, iir_a, x_mixed, [], 2);          % [N x T x T_trials]
    x_colored = x_colored / iir_std;

    noise = reshape(sigma_noise, N, 1, 1) .* x_colored;        % [N x T x T_trials]
end

% ------------------------------------------------------------------------- %

function fr = assemble_traces(params, cp, baseline, noise)
% ASSEMBLE_TRACES  Combine baseline, odor signal, and noise into full traces.
%
%   fr(i, t_frame, t_trial) = baseline(i)          for all t_frame
%                            + mu(i, t_trial)       for t_frame in odor window
%                            + noise(i, t_frame, t_trial)
%   then clipped to >= 0.

    N        = params.N;
    T        = params.T;
    T_trials = size(cp.mu, 2);
    f1       = params.odor_frames(1);
    f2       = params.odor_frames(2);

    % Initialize with per-neuron baseline: [N x T x T_trials]
    fr = repmat(baseline, 1, T, T_trials);

    % Add odor response during odor window only.
    % mu is [N x T_trials]; reshape to [N x 1 x T_trials] to broadcast over T.
    mu_t = reshape(cp.mu, N, 1, T_trials);
    fr(:, f1:f2, :) = fr(:, f1:f2, :) + mu_t;

    % Add noise and clip
    fr = max(fr + noise, 0);
end
