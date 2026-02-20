classdef GeometryBuilder
% GeometryBuilder  Constructs the geometric axes of the funnel model.
%
%   Builds the set of unit vectors {e_k, n_k, d} in R^N that define the
%   funnel manifold geometry, given a ToyParams object.
%
%   Identity axes e_k:  one per odor, drawn from a D-dimensional subspace U.
%                       Pairwise correlations <e_k, e_j> = rho_e exactly,
%                       set independently of D via a Gram matrix construction.
%
%   Novelty axes n_k:   n_k = alpha * e_k + sqrt(1-alpha^2) * w_k
%                       where w_k ⊥ span(U) with <w_k, w_j> = rho exactly.
%                       Since w_k ⊥ e_k, n_k is already a unit vector.
%                       Actual <n_k, n_j> = alpha^2 * <e_k, e_j> + (1-alpha^2) * rho.
%
%   Drift direction d:  single unit vector orthogonal to span(U).
%
%   Usage:
%     params = ToyParams('D', 20, 'rho', 0.3, 'rho_e', 0.2, 'alpha', 0.1);
%     geom   = GeometryBuilder(params);
%     geom.e          % [N x K] identity axes
%     geom.n          % [N x K] novelty axes
%     geom.d          % [N x 1] drift direction
%     geom.report()   % print diagnostic correlations
%
%   See also: ToyParams, CentralPatterns, ToyDataGenerator

    properties (SetAccess = private)

        e       double  % [N x K]  odor identity axes (unit vectors)
        n       double  % [N x K]  odor novelty axes  (unit vectors)
        d       double  % [N x 1]  shared drift direction (unit vector)
        U       double  % [N x D]  identity subspace basis (orthonormal columns)

        % Diagnostics: actual inner products after construction.
        % Useful for verifying that targets were met.
        e_corr  double  % [K x K]  <e_k, e_j>  (diagonal = 1 by construction)
        n_corr  double  % [K x K]  <n_k, n_j>  (diagonal = 1 by construction)
        en_corr double  % [K x 1]  <e_k, n_k>  (target: alpha for all k)

    end

    % --------------------------------------------------------------------- %
    methods

        function obj = GeometryBuilder(params)
        % GEOMETRYBUILDER  Construct geometric axes from a ToyParams object.
        %
        %   geom = GeometryBuilder(params)
        %
        %   Sets the RNG state from params.rng_seed before sampling.
        %   Call params.validate() first to check parameter consistency.

            params.validate();
            rng(params.rng_seed);

            N = params.N;
            K = params.K;
            D = params.D;

            % 1. Identity subspace: [N x D] orthonormal basis
            obj.U = build_identity_subspace(N, D, params.rotation_mix);

            % 2. Identity axes: [N x K] unit vectors in span(U) with <e_k,e_j> = rho_e
            obj.e = build_identity_axes(obj.U, K, params.rho_e);

            % 3. Novelty axes: [N x K] unit vectors with controlled correlations
            obj.n = build_novelty_axes(N, K, obj.U, obj.e, params.alpha, params.rho);

            % 4. Drift direction: [N x 1] unit vector orthogonal to span(U)
            obj.d = build_drift_direction(N, obj.U);

            % 5. Diagnostics
            obj.e_corr  = obj.e' * obj.e;
            obj.n_corr  = obj.n' * obj.n;
            obj.en_corr = diag(obj.e' * obj.n);
        end

        % -- Diagnostics ------------------------------------------------ %

        function report(obj)
        % REPORT  Print a summary of actual axis correlations after build.
            K = size(obj.e, 2);
            idx = logical(triu(ones(K), 1));  % upper-triangle mask, k~=j

            e_off  = obj.e_corr(idx);
            n_off  = obj.n_corr(idx);

            fprintf('GeometryBuilder diagnostics\n');
            fprintf('  Identity axes  <e_k, e_j>  (k~=j):  mean = %+.4f,  std = %.4f\n', ...
                    mean(e_off), std(e_off));
            fprintf('  Novelty axes   <n_k, n_j>  (k~=j):  mean = %+.4f,  std = %.4f\n', ...
                    mean(n_off), std(n_off));
            fprintf('  Within-odor    <e_k, n_k>:           mean = %+.4f,  std = %.4f\n', ...
                    mean(obj.en_corr), std(obj.en_corr));
            fprintf('  Drift          ||d|| = %.6f\n', norm(obj.d));
        end

    end

end

% ========================================================================= %
% Local helper functions
% ========================================================================= %

function U = build_identity_subspace(N, D, rotation_mix)
% BUILD_IDENTITY_SUBSPACE  Returns an [N x D] orthonormal basis for the
% identity subspace, interpolated between canonical and random orientation.
%
%   rotation_mix = 0: U = first D standard basis vectors (sparse neurons)
%   rotation_mix = 1: U = random orthonormal frame (dense mixed selectivity)
%   Intermediate: column-wise blend followed by QR re-orthonormalization.

    U_canonical = [eye(D); zeros(N - D, D)];        % [N x D]
    [U_random, ~] = qr(randn(N, D), 'econ');         % [N x D], Haar-distributed

    if rotation_mix == 0
        U = U_canonical;
    elseif rotation_mix == 1
        U = U_random;
    else
        U_blend = (1 - rotation_mix) * U_canonical + rotation_mix * U_random;
        [U, ~]  = qr(U_blend, 'econ');
    end
end

% ------------------------------------------------------------------------- %

function e = build_identity_axes(U, K, rho_e)
% BUILD_IDENTITY_AXES  Returns [N x K] unit vectors in span(U) with
%   <e_k, e_j> = rho_e exactly for all k ~= j (Gram matrix construction).
%
%   Steps:
%     1. Build Gram matrix G: G_{kj} = rho_e (k~=j), G_{kk} = 1.
%     2. Cholesky factor: G = L * L'  (regularized for boundary case).
%     3. Draw K orthonormal coordinate vectors in R^D via QR.
%     4. E_raw = U * V_e  → [N x K] orthonormal in R^N, within span(U).
%     5. e     = E_raw * L'  → <e_k, e_j> = G_{kj} = rho_e.
%
%   Proof: e'*e = L*(E_raw'*E_raw)*L' = L*I*L' = G.

    D = size(U, 2);

    % Gram matrix for identity axes
    G = rho_e * ones(K) + (1 - rho_e) * eye(K) + 1e-10 * eye(K);
    L = chol(G, 'lower');   % [K x K], G = L * L'

    % K orthonormal coordinate vectors in R^D
    [V_e, ~] = qr(randn(D, K), 'econ');   % [D x K], orthonormal
    V_e      = V_e(:, 1:K);

    % Map into neuron space via identity subspace
    E_raw = U * V_e;   % [N x K], orthonormal columns in R^N

    % Apply Cholesky factor to impose target correlations
    e = E_raw * L';    % [N x K], <e_k,e_j> = rho_e exactly
end

% ------------------------------------------------------------------------- %

function n = build_novelty_axes(N, K, U, e, alpha, rho)
% BUILD_NOVELTY_AXES  Returns [N x K] novelty axes with controlled correlations.
%
%   Construction:
%     w_k ⊥ span(U),  <w_k, w_j> = rho (exact, via Gram matrix Cholesky)
%     n_k = alpha * e_k + sqrt(1-alpha^2) * w_k
%
%   Since e_k ∈ span(U) and w_k ⊥ span(U):  <e_k, w_k> = 0
%   Therefore ||n_k|| = sqrt(alpha^2 + (1-alpha^2)) = 1  (no renorm needed).
%
%   Actual <n_k, n_j> = alpha^2 * <e_k, e_j> + (1-alpha^2) * rho.
%   Deviation from target rho grows with alpha^2 * mean|<e_k, e_j>|.

    % Gram matrix for w_k components: G_{kj} = rho (k~=j), G_{kk} = 1
    G   = rho * ones(K) + (1 - rho) * eye(K);

    % Regularize to handle rho at the PSD boundary (rho = -1/(K-1))
    G   = G + 1e-10 * eye(K);

    L   = chol(G, 'lower');  % [K x K], G = L * L'

    % K orthonormal vectors in the complement of span(U)
    V_raw     = randn(N, K);
    V_raw     = V_raw - U * (U' * V_raw);   % project out identity subspace
    [V, ~]    = qr(V_raw, 'econ');
    V         = V(:, 1:K);                  % [N x K], orthonormal, ⊥ U

    % W columns: w_k with <w_k, w_j> = G_{kj} = rho (k~=j), ||w_k|| = 1
    % Proof: W'*W = (V*L')'*(V*L') = L * (V'*V) * L' = L*I*L' = G
    W = V * L';  % [N x K]

    % Combine: n_k = alpha*e_k + sqrt(1-alpha^2)*w_k  (already unit norm)
    n = alpha * e + sqrt(1 - alpha^2) * W;
end

% ------------------------------------------------------------------------- %

function d = build_drift_direction(N, U)
% BUILD_DRIFT_DIRECTION  Returns an [N x 1] unit vector orthogonal to span(U).
%
%   The drift direction is in the complement of the identity subspace,
%   so the directed drift component is independent of odor identity axes.

    d_raw = randn(N, 1);
    d_raw = d_raw - U * (U' * d_raw);  % project out identity subspace
    d     = d_raw / norm(d_raw);
end
