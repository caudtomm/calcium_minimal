classdef CentralPatterns
% CentralPatterns  Computes the noiseless central population patterns.
%
%   Evaluates the funnel model's mean response vector for every trial t
%   (defined by stim_idx):
%
%     mu(:,t) = A*(1+lambda_A*(1-exp(-(r-1)/tau_A))) * e_k(r)   [identity]
%             + B*exp(-(r-1)/tau)                    * n_k       [novelty]
%             + gamma*(r-1)                          * d_k       [drift]
%             + C*[cos(theta*pi/2*(t-1))*u_s
%                + sin(theta*pi/2*(t-1))*v_s]                   [baseline]
%
%   where k = stim_idx(t) and r = number of times odor k has appeared
%   in stim_idx(1..t) (1-indexed repetition count).
%
%   e_k(r): identity axes, may vary across reps when eta > 0.
%   n_k, d_k: fixed (built from rep-1 geometry).
%   u_s, v_s: baseline axes, orthogonal to all others.
%   Novelty uses tau; identity growth (lambda_A, eta) uses tau_A.
%   ||s_t|| = C exactly.  theta=0: constant offset; theta=1: consecutive
%   trials orthogonal.
%
%   mu is the expected population firing rate (Hz) during the odor window,
%   relative to zero.  The NoiseModel adds a per-neuron baseline so that
%   total firing rates remain positive.
%
%   Null model (B=0, gamma=0, lambda_A=0, C=0): mu(:,t) identical for all
%   trials presenting the same odor at the same repetition.
%
%   Usage:
%     params = ToyParams();
%     geom   = GeometryBuilder(params);
%     cp     = CentralPatterns(params, geom);
%     cp.mu           % [N x T_trials] central patterns
%     cp.report()     % print summary statistics
%
%   See also: ToyParams, GeometryBuilder, NoiseModel, ToyDataGenerator

    properties (SetAccess = private)

        mu               double  % [N x T_trials]  central patterns (Hz)
        stim_idx         double  % [1 x T_trials]  trial presentation order

        % Per-repetition scalar weights (indexed by rep number, length R_max)
        identity_weights double  % [R_max x 1]  A*(1+lambda_A*(1-exp(-(r-1)/tau_A)))
        novelty_weights  double  % [R_max x 1]  B*exp(-(r-1)/tau)
        drift_offsets    double  % [R_max x 1]  gamma*(r-1)
        rho_e_schedule   double  % [R_max x 1]  rho_e+eta*(1-exp(-(r-1)/tau_A))

    end

    % --------------------------------------------------------------------- %
    methods

        function obj = CentralPatterns(params, geom)
        % CENTRALPATTERNS  Build central patterns from params and geometry.
        %
        %   cp = CentralPatterns(params, geom)
        %
        %   params : ToyParams
        %   geom   : GeometryBuilder

            N = params.N;
            K = params.K;

            obj.stim_idx = params.get_stim_idx();   % [1 x T_trials]
            T_trials     = numel(obj.stim_idx);

            % Max repetition count across all odors
            rep_hist = histcounts(obj.stim_idx, 0.5 : K + 0.5);  % [1 x K]
            R_max    = max(rep_hist);

            % Pre-compute per-rep weights up to R_max
            r_idx         = (0 : R_max-1)';
            exp_decay_nov = exp(-r_idx / params.tau);
            exp_decay_id  = exp(-r_idx / params.tau_A);

            obj.identity_weights = params.A * (1 + params.lambda_A * (1 - exp_decay_id));
            obj.novelty_weights  = params.B * exp_decay_nov;
            obj.drift_offsets    = params.gamma * r_idx;
            obj.rho_e_schedule   = params.rho_e + params.eta * (1 - exp_decay_id);

            % Allocate output: [N x T_trials]
            obj.mu    = zeros(N, T_trials);

            % Angular drift rate for baseline (radians per trial step)
            omega     = params.theta * (pi / 2);

            rep_count = zeros(1, K);  % how many times each odor has appeared

            for t = 1:T_trials
                k = obj.stim_idx(t);
                rep_count(k) = rep_count(k) + 1;
                r = rep_count(k);

                iw = obj.identity_weights(r);
                nw = obj.novelty_weights(r);
                dw = obj.drift_offsets(r);

                % Per-rep identity axes: fixed E_raw, rep-dependent Gram mixing
                rho_r = obj.rho_e_schedule(r);
                G_r   = rho_r * ones(K) + (1 - rho_r) * eye(K) + 1e-10 * eye(K);
                L_r   = chol(G_r, 'lower');
                e_r   = geom.E_raw * L_r';  % [N x K]

                e_k = e_r(:, k);        % [N x 1]
                n_k = geom.n(:, k);     % [N x 1]
                d_k = geom.d(:, k);     % [N x 1]

                % Baseline: rotates in (u_s, v_s) plane at angular speed omega
                phase = omega * (t - 1);
                s_t   = params.C * (cos(phase) * geom.u_s + sin(phase) * geom.v_s);

                obj.mu(:, t) = iw * e_k + nw * n_k + dw * d_k + s_t;
            end
        end

        % -- Diagnostics ------------------------------------------------ %

        function report(obj)
        % REPORT  Print a summary of central pattern statistics.
            [~, T_trials] = size(obj.mu);
            R_max = numel(obj.novelty_weights);

            fprintf('CentralPatterns diagnostics\n');
            fprintf('  Trial sequence: T_trials = %d\n', T_trials);

            fprintf('  rho_e schedule  (rho_e + eta*(1-exp(-(r-1)/tau_A))):\n');
            for r = 1:R_max
                fprintf('    rep %d:  %.4f\n', r, obj.rho_e_schedule(r));
            end
            fprintf('  Identity weights  A*(1+lambda_A*(1-exp(-(r-1)/tau_A))):\n');
            for r = 1:R_max
                fprintf('    rep %d:  %.4f Hz\n', r, obj.identity_weights(r));
            end
            fprintf('  Novelty weights  B*exp(-(r-1)/tau):\n');
            for r = 1:R_max
                fprintf('    rep %d:  %.4f Hz\n', r, obj.novelty_weights(r));
            end
            fprintf('  Drift offsets  gamma*(r-1):\n');
            for r = 1:R_max
                fprintf('    rep %d:  %.4f Hz\n', r, obj.drift_offsets(r));
            end

            fprintf('  Per-trial mean pattern norm  ||mu(:,t)||:\n');
            norms = vecnorm(obj.mu, 2, 1);   % [1 x T_trials]
            fprintf('    mean = %.4f,  std = %.4f,  min = %.4f,  max = %.4f\n', ...
                    mean(norms), std(norms), min(norms), max(norms));

            fprintf('  Fraction of (neuron, trial) entries < 0:  %.2f%%\n', ...
                    100 * mean(obj.mu(:) < 0));
        end

    end

end
