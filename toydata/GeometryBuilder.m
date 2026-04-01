classdef GeometryBuilder
% GeometryBuilder  Constructs the geometric axes of the funnel model.
%
%   Builds the set of unit vectors {e_k, n_k, d} in R^N that define the
%   funnel manifold geometry, given a ToyParams object.
%
%   Identity axes e_k:  one per odor, drawn from a D-dimensional subspace U.
%                       Pairwise correlations <e_k, e_j> = rho_e exactly at
%                       rep 1, set via a Gram matrix construction.
%                       E_raw (the underlying orthonormal frame) is stored so
%                       that CentralPatterns can recompute e_k per repetition
%                       when eta > 0, by applying a new Cholesky factor to the
%                       same E_raw without any N-D rotation.
%
%   Novelty axes n_k:   n_k = alpha * e_k + sqrt(1-alpha^2) * w_k
%                       where w_k ⊥ span(U) with <w_k, w_j> = rho exactly.
%                       Since w_k ⊥ e_k, n_k is already a unit vector.
%                       Actual <n_k, n_j> = alpha^2 * <e_k, e_j> + (1-alpha^2) * rho.
%
%   Drift directions d_k: one per odor, orthogonal to span(U).
%                         Pairwise correlations <d_k, d_j> = rho_d exactly.
%
%   Usage:
%     params = ToyParams('D', 20, 'rho', 0.3, 'rho_e', 0.2, 'rho_d', 0.5, 'alpha', 0.1);
%     geom   = GeometryBuilder(params);
%     geom.e          % [N x K] identity axes
%     geom.n          % [N x K] novelty axes
%     geom.d          % [N x 1] drift direction
%     geom.report()   % print diagnostic correlations
%
%   See also: ToyParams, CentralPatterns, ToyDataGenerator

    properties (SetAccess = private)

        e       double  % [N x K]  odor identity axes at rep 1 (unit vectors)
        E_raw   double  % [N x K]  orthonormal frame underlying e (U * V_e);
                        %          fixed across reps; used by CentralPatterns
                        %          to recompute e per repetition when eta > 0.
        n       double  % [N x K]  odor novelty axes  (unit vectors)
        d       double  % [N x K]  odor-specific drift directions (unit vectors)
        u_s     double  % [N x 1]  baseline axis (unit vector, ⊥ all others)
        v_s     double  % [N x 1]  baseline rotation axis (unit vector, ⊥ all + u_s)
        U       double  % [N x D]  identity subspace basis (orthonormal columns)

        % Diagnostics: actual inner products after construction.
        % Useful for verifying that targets were met.
        e_corr  double  % [K x K]  <e_k, e_j>  (diagonal = 1 by construction)
        n_corr  double  % [K x K]  <n_k, n_j>  (diagonal = 1 by construction)
        d_corr  double  % [K x K]  <d_k, d_j>  (diagonal = 1 by construction)
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
            [obj.e, obj.E_raw] = build_identity_axes(obj.U, K, params.rho_e);

            % 3. Novelty axes: [N x K] unit vectors with controlled correlations
            obj.n = build_novelty_axes(N, K, obj.U, obj.e, params.alpha, params.rho);

            % 4. Drift directions: [N x K] unit vectors orthogonal to span(U)
            obj.d = build_drift_directions(N, K, obj.U, params.rho_d);

            % 5. Baseline axes: unit vectors orthogonal to span(U, n, d)
            all_axes = [obj.U, obj.n, obj.d];
            obj.u_s  = build_baseline_axis(N, all_axes);
            obj.v_s  = build_baseline_axis(N, [all_axes, obj.u_s]);

            % 6. Diagnostics
            obj.e_corr  = obj.e' * obj.e;
            obj.n_corr  = obj.n' * obj.n;
            obj.d_corr  = obj.d' * obj.d;
            obj.en_corr = diag(obj.e' * obj.n);
        end

        % -- Diagnostics ------------------------------------------------ %

        function report(obj)
        % REPORT  Print a summary of actual axis correlations after build.
            K = size(obj.e, 2);
            idx = logical(triu(ones(K), 1));  % upper-triangle mask, k~=j

            e_off  = obj.e_corr(idx);
            n_off  = obj.n_corr(idx);
            d_off  = obj.d_corr(idx);

            fprintf('GeometryBuilder diagnostics\n');
            fprintf('  Identity axes  <e_k, e_j>  (k~=j):  mean = %+.4f,  std = %.4f\n', ...
                    mean(e_off), std(e_off));
            fprintf('  Novelty axes   <n_k, n_j>  (k~=j):  mean = %+.4f,  std = %.4f\n', ...
                    mean(n_off), std(n_off));
            fprintf('  Drift dirs     <d_k, d_j>  (k~=j):  mean = %+.4f,  std = %.4f\n', ...
                    mean(d_off), std(d_off));
            fprintf('  Within-odor    <e_k, n_k>:           mean = %+.4f,  std = %.4f\n', ...
                    mean(obj.en_corr), std(obj.en_corr));
            fprintf('  Baseline axes  ||u_s|| = %.6f,  ||v_s|| = %.6f\n', ...
                    norm(obj.u_s), norm(obj.v_s));
            fprintf('                 <u_s, v_s> = %+.2e  (target: 0)\n', ...
                    dot(obj.u_s, obj.v_s));
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

function [e, E_raw] = build_identity_axes(U, K, rho_e)
% BUILD_IDENTITY_AXES  Returns [N x K] identity axes with <e_k,e_j> = rho_e,
%   plus the underlying orthonormal frame E_raw for per-rep updates.
%
%   Steps:
%     1. Build Gram matrix G: G_{kj} = rho_e (k~=j), G_{kk} = 1.
%     2. Cholesky factor: G = L * L'  (regularized for boundary case).
%     3. Draw K orthonormal coordinate vectors in R^D via QR.
%     4. E_raw = U * V_e  → [N x K] orthonormal in R^N, within span(U).
%     5. e     = E_raw * L'  → <e_k, e_j> = G_{kj} = rho_e.
%
%   Proof: e'*e = L*(E_raw'*E_raw)*L' = L*I*L' = G.
%
%   E_raw is returned so that CentralPatterns can recompute e_k for a
%   different rho_e_eff(r) without any N-dimensional rotation:
%     e_r = E_raw * chol(G(rho_e_eff(r)), 'lower')'

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

function d = build_drift_directions(N, K, U, rho_d)
% BUILD_DRIFT_DIRECTIONS  Returns [N x K] odor-specific drift directions
%   orthogonal to span(U), with <d_k, d_j> = rho_d exactly (k ~= j).
%
%   Construction mirrors build_novelty_axes (without the alpha mixing):
%     V: K orthonormal vectors in complement of span(U)
%     G = rho_d*ones(K) + (1-rho_d)*eye(K)  → Gram matrix
%     d = V * chol(G, 'lower')'
%   Proof: d'*d = L*(V'V)*L' = L*I*L' = G.

    G   = rho_d * ones(K) + (1 - rho_d) * eye(K) + 1e-10 * eye(K);
    L   = chol(G, 'lower');

    V_raw = randn(N, K);
    V_raw = V_raw - U * (U' * V_raw);  % project out identity subspace
    [V, ~] = qr(V_raw, 'econ');
    V      = V(:, 1:K);                % [N x K], orthonormal, ⊥ U

    d = V * L';   % [N x K], <d_k, d_j> = rho_d exactly
end

% ------------------------------------------------------------------------- %

function u = build_baseline_axis(N, all_axes)
% BUILD_BASELINE_AXIS  Returns a unit vector orthogonal to col(all_axes).
%
%   Computes the QR orthonormal basis of all_axes, then projects a random
%   vector onto its complement and normalises.  Exact orthogonality (not
%   approximate) for u_s and v_s relative to span(U, n, d).

    [Q, ~]  = qr(all_axes, 'econ');   % orthonormal basis for col(all_axes)
    u_raw   = randn(N, 1);
    u_raw   = u_raw - Q * (Q' * u_raw);
    u       = u_raw / norm(u_raw);
end
