classdef ToyDataGenerator
% ToyDataGenerator  Top-level orchestrator for the funnel generative model.
%
%   Runs the full pipeline in sequence:
%     ToyParams -> GeometryBuilder -> CentralPatterns -> NoiseModel
%
%   and packages the result for saving and downstream analysis.
%
%   Usage (single dataset):
%     params = ToyParams('tau', 2, 'rho', 0.3);
%     tg     = ToyDataGenerator(params);
%     tg.saveToFile('/path/to/outdir');
%     tg.report();
%
%   Usage (batch / HPC):
%     ToyDataGenerator.generate_batch(params, outdir, n_replicas);
%     % or, for a parameter sweep:
%     ToyDataGenerator.generate_batch(params_cell, outdir);
%
%   Output .mat file contains:
%     firing_rates   [N x T x K x R]  double   firing rates (Hz)
%     params         struct                     all ToyParams fields
%     metadata       struct                     timing, labels, model type
%     geometry       struct                     axes e, n, d, subspace U
%     patterns       struct                     mu, novelty_weights, drift_offsets
%
%   See also: ToyParams, GeometryBuilder, CentralPatterns, NoiseModel, ToyDataLoader

    properties (SetAccess = private)
        params  ToyParams        % input parameters
        geom    GeometryBuilder  % geometric axes
        cp      CentralPatterns  % noiseless central patterns
        nm      NoiseModel       % noise model + final firing rates
    end

    % --------------------------------------------------------------------- %
    methods

        function obj = ToyDataGenerator(params)
        % TOYDATAGENERATOR  Run the full generative pipeline.
        %
        %   tg = ToyDataGenerator(params)
        %
        %   Sets the global RNG state from params.rng_seed before generation.

            params.validate();
            rng(params.rng_seed);

            obj.params = params;
            obj.geom   = GeometryBuilder(params);
            obj.cp     = CentralPatterns(params, obj.geom);
            obj.nm     = NoiseModel(params, obj.cp);
        end

        % -- I/O -------------------------------------------------------- %

        function filepath = saveToFile(obj, outdir, filename)
        % SAVETOFILE  Save generated data to a .mat file (v7.3 / HDF5).
        %
        %   filepath = tg.saveToFile(outdir)
        %   filepath = tg.saveToFile(outdir, filename)
        %
        %   If filename is omitted, a descriptive name is generated from
        %   key parameter values (see default_filename).

            if nargin < 2 || isempty(outdir)
                outdir = pwd;
            end
            if nargin < 3 || isempty(filename)
                filename = obj.default_filename();
            end

            if ~exist(outdir, 'dir')
                mkdir(outdir);
            end

            filepath = fullfiletol(outdir, filename);
            s = obj.toSaveStruct();

            % Use v7.3 (HDF5) for large arrays and Python/mat73 compatibility
            save(filepath, '-struct', 's', '-v7.3');
            fprintf('Saved: %s\n', filepath);
        end

        function s = toSaveStruct(obj)
        % TOSAVESTRUCT  Package all outputs into a plain struct for saving.

            s.params   = obj.params.toStruct();
            s.metadata = build_metadata(obj);

            % Primary output
            s.firing_rates = obj.nm.firing_rates;   % [N x T x T_trials]

            % Geometric axes (small; useful for validation and the loader)
            s.geometry.e   = obj.geom.e;            % [N x K]
            s.geometry.n   = obj.geom.n;            % [N x K]
            s.geometry.d   = obj.geom.d;            % [N x K]
            s.geometry.U   = obj.geom.U;            % [N x D]
            s.geometry.u_s = obj.geom.u_s;          % [N x 1]
            s.geometry.v_s = obj.geom.v_s;          % [N x 1]

            % Noiseless central patterns (small)
            s.patterns.mu               = obj.cp.mu;                % [N x T_trials]
            s.patterns.stim_idx         = obj.cp.stim_idx;          % [1 x T_trials]
            s.patterns.identity_weights = obj.cp.identity_weights;  % [R_max x 1]
            s.patterns.novelty_weights  = obj.cp.novelty_weights;   % [R_max x 1]
            s.patterns.drift_offsets    = obj.cp.drift_offsets;     % [R_max x 1]
            s.patterns.rho_e_schedule   = obj.cp.rho_e_schedule;    % [R_max x 1]

            % Noise diagnostics
            s.noise.baseline    = obj.nm.baseline;     % [N x 1]
            s.noise.sigma_noise = obj.nm.sigma_noise;  % [N x 1]
            s.noise.ac_factor   = obj.nm.ac_factor;    % scalar
        end

        % -- Diagnostics ------------------------------------------------ %

        function report(obj)
        % REPORT  Print a full diagnostic summary of the generated dataset.
            fprintf('=== ToyDataGenerator report ===\n');
            fprintf('Model type: %s\n\n', obj.nm_model_type());
            obj.geom.report();
            fprintf('\n');
            obj.cp.report();
            fprintf('\n');
            obj.nm.report(obj.params);
        end

    end

    % --------------------------------------------------------------------- %
    methods (Access = private)

        function fname = default_filename(obj)
        % DEFAULT_FILENAME  Generate a descriptive filename from key params.
            p = obj.params;
            fname = sprintf('toydata_%s_N%d_K%d_R%d_D%d_tau%.1f_rho%.2f_seed%d.mat', ...
                obj.nm_model_type(), p.N, p.K, p.R, p.D, p.tau, p.rho, p.rng_seed);
        end

        function s = nm_model_type(obj)
        % NM_MODEL_TYPE  Returns 'null' or 'hypothesis' string.
            if obj.params.B == 0 && obj.params.gamma == 0
                s = 'null';
            else
                s = 'hypothesis';
            end
        end

    end

    % --------------------------------------------------------------------- %
    methods (Static)

        function generate_batch(params, outdir, n_replicas)
        % GENERATE_BATCH  Generate multiple datasets, saving each to outdir.
        %
        %   ToyDataGenerator.generate_batch(params, outdir, n_replicas)
        %     Generates n_replicas datasets from a single ToyParams, using
        %     seeds rng_seed, rng_seed+1, ..., rng_seed+n_replicas-1.
        %
        %   ToyDataGenerator.generate_batch(params_cell, outdir)
        %     Generates one dataset per entry in a cell array of ToyParams.
        %     Useful for parameter sweeps.

            if iscell(params)
                param_list  = params;
                n_replicas  = numel(param_list);
            else
                if nargin < 3
                    error('ToyDataGenerator.generate_batch: n_replicas required when params is a ToyParams.');
                end
                base_seed  = params.rng_seed;
                param_list = cell(n_replicas, 1);
                for i = 1:n_replicas
                    p           = params;
                    p.rng_seed  = base_seed + i - 1;
                    param_list{i} = p;
                end
            end

            if ~exist(outdir, 'dir')
                mkdir(outdir);
            end

            for i = 1:n_replicas
                fprintf('[%d / %d]  seed = %d\n', i, n_replicas, param_list{i}.rng_seed);
                tg = ToyDataGenerator(param_list{i});
                tg.saveToFile(outdir);
            end

            fprintf('Done. %d file(s) written to: %s\n', n_replicas, outdir);
        end

    end

end

% ========================================================================= %
% Local helper functions
% ========================================================================= %

function meta = build_metadata(obj)
% BUILD_METADATA  Assemble a self-contained metadata struct for the loader.
%
%   Contains everything ToyDataLoader needs without requiring ToyParams.

    p = obj.params;

    meta.N              = p.N;
    meta.K              = p.K;
    meta.R              = p.R;
    meta.T              = p.T;
    meta.fs             = p.fs;
    meta.time           = (0 : p.T - 1) / p.fs;       % [1 x T] seconds
    meta.t_odor_start   = p.t_odor_start;
    meta.t_odor_end     = p.t_odor_end;
    meta.odor_frames    = p.odor_frames;               % [start, end] 1-indexed
    meta.stimulus_names = p.get_stimulus_names();      % {1 x K} cell of char
    meta.stim_idx       = obj.cp.stim_idx;             % [1 x T_trials]
    meta.T_trials       = numel(obj.cp.stim_idx);

    if p.B == 0 && p.gamma == 0 && p.lambda_A == 0 && p.C == 0
        meta.model_type = 'null';
    else
        meta.model_type = 'hypothesis';
    end

    meta.baseline       = obj.nm.baseline;   % [N x 1] per-neuron baseline (Hz)
end
