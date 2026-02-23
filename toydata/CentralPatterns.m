classdef CentralPatterns
% CentralPatterns  Computes the noiseless central population patterns.
%
%   Evaluates the funnel model's mean response vector for every (odor, rep)
%   combination:
%
%     mu(:,k,r) = A*(1 + lambda_A*(1-exp(-(r-1)/tau))) * e_k
%               + B*exp(-(r-1)/tau)                    * n_k
%               + gamma*(r-1)                          * d
%               |________identity (growing)___________|
%                         |_____novelty (decaying)_____|  |__drift__|
%
%   Identity and novelty share the time constant tau: as novelty decays,
%   identity grows toward A*(1+lambda_A).  lambda_A=0 recovers the original
%   flat identity model.
%
%   mu is the expected population firing rate (Hz) during the odor window,
%   relative to zero.  It may be negative for some neurons; the NoiseModel
%   adds a per-neuron baseline so that total firing rates remain positive.
%
%   Null model (B=0, gamma=0, lambda_A=0): mu(:,:,r) is identical for all r.
%
%   Usage:
%     params = ToyParams();
%     geom   = GeometryBuilder(params);
%     cp     = CentralPatterns(params, geom);
%     cp.mu           % [N x K x R] central patterns
%     cp.report()     % print summary statistics
%
%   See also: ToyParams, GeometryBuilder, NoiseModel, ToyDataGenerator

    properties (SetAccess = private)

        mu              double  % [N x K x R]  central patterns (Hz)

        % Per-repetition scalar weights (useful for plotting trajectories)
        identity_weights double % [R x 1]  A*(1+lambda_A*(1-exp(-(r-1)/tau))) for r=1..R
        novelty_weights  double % [R x 1]  B*exp(-(r-1)/tau)                  for r=1..R
        drift_offsets    double % [R x 1]  gamma*(r-1)                         for r=1..R
        rho_e_schedule   double % [R x 1]  rho_e+eta*(1-exp(-(r-1)/tau))      for r=1..R

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
            R = params.R;

            r_idx = (0 : R-1)';                        % [R x 1], 0-indexed
            exp_decay = exp(-r_idx / params.tau);      % [R x 1], shared decay envelope

            obj.identity_weights = params.A * (1 + params.lambda_A * (1 - exp_decay));
            obj.novelty_weights  = params.B * exp_decay;
            obj.drift_offsets    = params.gamma * r_idx;
            obj.rho_e_schedule   = params.rho_e + params.eta * (1 - exp_decay);

            % Allocate output: [N x K x R]
            obj.mu = zeros(N, K, R);

            for r = 1:R
                iw = obj.identity_weights(r);  % scalar
                nw = obj.novelty_weights(r);   % scalar
                dw = obj.drift_offsets(r);     % scalar

                % Per-rep identity axes: same orthonormal frame E_raw, new
                % Gram mixing for rho_e_eff(r).  n_k is kept fixed (rep-1).
                rho_r = obj.rho_e_schedule(r);
                G_r   = rho_r * ones(K) + (1 - rho_r) * eye(K) + 1e-10 * eye(K);
                L_r   = chol(G_r, 'lower');
                e_r   = geom.E_raw * L_r';     % [N x K]

                % iw*[N x K] + nw*[N x K] + dw*[N x 1]  (column broadcast)
                obj.mu(:, :, r) = iw * e_r  + nw * geom.n  + dw * geom.d;
            end
        end

        % -- Diagnostics ------------------------------------------------ %

        function report(obj)
        % REPORT  Print a summary of central pattern statistics.
            [~, K, R] = size(obj.mu);

            fprintf('CentralPatterns diagnostics\n');
            fprintf('  rho_e schedule  (rho_e + eta*(1-exp(-(r-1)/tau))):\n');
            for r = 1:R
                fprintf('    rep %d:  %.4f\n', r, obj.rho_e_schedule(r));
            end
            fprintf('  Identity weights  A*(1+lambda_A*(1-exp(-(r-1)/tau))):\n');
            for r = 1:R
                fprintf('    rep %d:  %.4f Hz\n', r, obj.identity_weights(r));
            end
            fprintf('  Novelty weights  B*exp(-(r-1)/tau):\n');
            for r = 1:R
                fprintf('    rep %d:  %.4f Hz\n', r, obj.novelty_weights(r));
            end
            fprintf('  Drift offsets  gamma*(r-1):\n');
            for r = 1:R
                fprintf('    rep %d:  %.4f Hz\n', r, obj.drift_offsets(r));
            end

            fprintf('  Per-rep mean pattern norm  ||mu(:,k,r)||  (averaged over k):\n');
            for r = 1:R
                norms = vecnorm(obj.mu(:, :, r), 2, 1);  % [1 x K]
                fprintf('    rep %d:  mean = %.4f,  std = %.4f\n', ...
                        r, mean(norms), std(norms));
            end

            fprintf('  Fraction of (neuron, odor, rep) entries < 0:  %.2f%%\n', ...
                    100 * mean(obj.mu(:) < 0));
        end

    end

end
