classdef CentralPatterns
% CentralPatterns  Computes the noiseless central population patterns.
%
%   Evaluates the funnel model's mean response vector for every (odor, rep)
%   combination:
%
%     mu(:, k, r) = A * e_k  +  B * exp(-(r-1)/tau) * n_k  +  gamma*(r-1) * d
%                  |_identity_|  |_______novelty decay_______|  |___drift___|
%
%   mu is the expected population firing rate (Hz) during the odor window,
%   relative to zero.  It may be negative for some neurons; the NoiseModel
%   adds a per-neuron baseline so that total firing rates remain positive.
%
%   Null model (B=0, gamma=0): mu(:,:,r) is identical for all r.
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
        novelty_weights double  % [R x 1]  B * exp(-(r-1)/tau) for r = 1..R
        drift_offsets   double  % [R x 1]  gamma * (r-1)        for r = 1..R

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

            obj.novelty_weights = params.B     * exp(-r_idx / params.tau);
            obj.drift_offsets   = params.gamma * r_idx;

            % Allocate output: [N x K x R]
            obj.mu = zeros(N, K, R);

            % Identity component is the same for all reps: [N x K]
            identity_part = params.A * geom.e;

            for r = 1:R
                nw = obj.novelty_weights(r);   % scalar
                dw = obj.drift_offsets(r);     % scalar

                % [N x K] + nw*[N x K] + dw*[N x 1]  (column broadcast)
                obj.mu(:, :, r) = identity_part  + nw * geom.n  + dw * geom.d;
            end
        end

        % -- Diagnostics ------------------------------------------------ %

        function report(obj)
        % REPORT  Print a summary of central pattern statistics.
            [~, K, R] = size(obj.mu);

            fprintf('CentralPatterns diagnostics\n');
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
