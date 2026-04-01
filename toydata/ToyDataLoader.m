classdef ToyDataLoader
% ToyDataLoader  Adapter from ToyDataGenerator .mat files to the analysis pipeline.
%
%   Reads .mat files produced by ToyDataGenerator and returns either
%   individual ToyActivityTraces objects or a fully assembled Experiment
%   that can be passed directly to ExperimentViewer, DataFilter,
%   ModeSelector, and GCMC_Analysis.
%
%   Each .mat file corresponds to one "subject" (one replica of the model).
%   Multiple files with different rng_seeds map to multiple subjects in one
%   experiment — the intended usage for statistical comparison with bio data.
%
%   Usage:
%     % Single file → ToyActivityTraces
%     at = ToyDataLoader.loadTraces('toydata_hypothesis_...mat', 'trained');
%
%     % Multiple files → Experiment
%     exp = ToyDataLoader.loadExperiment(filepaths, 'trained', 'my_toy_exp');
%
%     % Whole directory → Experiment
%     exp = ToyDataLoader.fromDir('/path/to/toydata/', 'trained');
%
%     % Then use exactly like biological data:
%     viewer = ExperimentViewer(exp);
%     viewer.dataFilter.subjectGroup = 'trained';
%     GCMC_Analysis(viewer).outputDataFiles('odor', outdir);
%
%   Note on Experiment construction:
%     Experiment requires a valid spreadsheet path ({isfile} validator).
%     ToyDataLoader writes a minimal temporary CSV, constructs the object,
%     then deletes the file and sets record.tabpath to a meaningful path.
%     Downstream code that calls fileparts(experiment.record.tabpath) for
%     naming will receive the exp_name, which is correct.
%
%   See also: ToyActivityTraces, ToyDataGenerator, ExperimentViewer

    % --------------------------------------------------------------------- %
    methods (Static)

        function at = loadTraces(matfile, group)
        % LOADTRACES  Load a single .mat file and return a ToyActivityTraces.
        %
        %   at = ToyDataLoader.loadTraces(matfile)
        %   at = ToyDataLoader.loadTraces(matfile, group)
        %
        %   matfile : char  path to a ToyDataGenerator .mat file
        %   group   : char  group label (default: inferred from params or 'toy')

            s = load(matfile);

            if nargin < 2 || isempty(group)
                group = infer_group(s);
            end

            at = ToyActivityTraces(s, group);
        end

        % ----------------------------------------------------------------- %

        function exp = loadExperiment(filepaths, group, exp_name)
        % LOADEXPERIMENT  Assemble an Experiment from a list of .mat files.
        %
        %   exp = ToyDataLoader.loadExperiment(filepaths, group)
        %   exp = ToyDataLoader.loadExperiment(filepaths, group, exp_name)
        %
        %   filepaths : char or {1xM cell} of .mat file paths
        %   group     : char  group label applied to all subjects
        %   exp_name  : char  used for experiment naming and output paths
        %                     (default: 'toy_experiment')
        %
        %   Returns an Experiment whose traces{i} are ToyActivityTraces.

            if ~iscell(filepaths)
                filepaths = {filepaths};
            end
            if nargin < 2 || isempty(group)
                group = 'toy';
            end
            if nargin < 3 || isempty(exp_name)
                exp_name = 'toy_experiment';
            end

            n = numel(filepaths);

            % Build subject table: one row per file
            names = arrayfun(@(i) sprintf('toy%03d', i), 1:n, 'UniformOutput', false);
            subjectTab = table(names', repmat({group}, n, 1), ...
                               'VariableNames', {'name', 'group'});

            % Experiment requires an existing file — write a minimal temp CSV
            tempfile = [tempname(), '.csv'];
            writetable(subjectTab, tempfile);

            try
                exp = Experiment(tempfile, 'toydata', '');
            catch err
                delete(tempfile);
                rethrow(err);
            end
            delete(tempfile);

            % Override tabpath with a meaningful name (used by GCMC_Analysis
            % to derive output directory names via fileparts).
            outdir_ref = fileparts(filepaths{1});
            exp.record.tabpath = fullfiletol(outdir_ref, [exp_name, '.csv']);
            exp.name       = string(exp_name);
            exp.subjectTab = subjectTab;

            % Load each file into a ToyActivityTraces
            fprintf('ToyDataLoader: loading %d file(s)...\n', n);
            for i = 1:n
                exp.traces{i} = ToyDataLoader.loadTraces(filepaths{i}, group);
                fprintf('  [%d/%d]  %s\n', i, n, filepaths{i});
            end
            fprintf('Done.\n');
        end

        % ----------------------------------------------------------------- %

        function exp = fromDir(dirpath, group, pattern, exp_name)
        % FROMDIR  Load all matching .mat files from a directory.
        %
        %   exp = ToyDataLoader.fromDir(dirpath)
        %   exp = ToyDataLoader.fromDir(dirpath, group)
        %   exp = ToyDataLoader.fromDir(dirpath, group, pattern)
        %   exp = ToyDataLoader.fromDir(dirpath, group, pattern, exp_name)
        %
        %   dirpath  : char  directory containing ToyDataGenerator .mat files
        %   group    : char  group label (default: 'toy')
        %   pattern  : char  glob pattern (default: 'toydata_*.mat')
        %   exp_name : char  experiment name (default: folder name of dirpath)

            if nargin < 2 || isempty(group);    group   = 'toy';            end
            if nargin < 3 || isempty(pattern);  pattern = 'toydata_*.mat'; end
            if nargin < 4 || isempty(exp_name)
                [~, exp_name] = fileparts(dirpath);
                if isempty(exp_name); exp_name = 'toy_experiment'; end
            end

            files = dir(fullfiletol(dirpath, pattern));
            if isempty(files)
                error('ToyDataLoader.fromDir: no files matching "%s" in:\n  %s', ...
                      pattern, dirpath);
            end

            filepaths = arrayfun(@(f) fullfiletol(f.folder, f.name), ...
                                 files, 'UniformOutput', false);

            exp = ToyDataLoader.loadExperiment(filepaths, group, exp_name);
        end

        % ----------------------------------------------------------------- %

        function exp = loadNullAndHypothesis(null_dir, hyp_dir, null_group, hyp_group)
        % LOADNULLANDHYPOTHESIS  Combine null and hypothesis replicas in one Experiment.
        %
        %   exp = ToyDataLoader.loadNullAndHypothesis(null_dir, hyp_dir)
        %   exp = ToyDataLoader.loadNullAndHypothesis(null_dir, hyp_dir, null_group, hyp_group)
        %
        %   Loads all 'toydata_null_*.mat' from null_dir and all
        %   'toydata_hypothesis_*.mat' from hyp_dir into a single Experiment,
        %   with subjects labelled by their respective groups.
        %   Useful for comparing null vs hypothesis in one ExperimentViewer.

            if nargin < 3 || isempty(null_group); null_group = 'null';       end
            if nargin < 4 || isempty(hyp_group);  hyp_group  = 'hypothesis'; end

            null_files = dir(fullfiletol(null_dir, 'toydata_null_*.mat'));
            hyp_files  = dir(fullfiletol(hyp_dir,  'toydata_hypothesis_*.mat'));

            if isempty(null_files)
                error('ToyDataLoader: no null .mat files found in: %s', null_dir);
            end
            if isempty(hyp_files)
                error('ToyDataLoader: no hypothesis .mat files found in: %s', hyp_dir);
            end

            null_paths = arrayfun(@(f) fullfiletol(f.folder, f.name), ...
                                  null_files, 'UniformOutput', false);
            hyp_paths  = arrayfun(@(f) fullfiletol(f.folder, f.name), ...
                                  hyp_files,  'UniformOutput', false);

            all_paths = [null_paths(:); hyp_paths(:)];
            all_groups = [repmat({null_group}, numel(null_paths), 1); ...
                          repmat({hyp_group},  numel(hyp_paths),  1)];

            % Build combined subject table
            n = numel(all_paths);
            names = arrayfun(@(i) sprintf('toy%03d', i), 1:n, 'UniformOutput', false);
            subjectTab = table(names', all_groups, 'VariableNames', {'name', 'group'});

            tempfile = [tempname(), '.csv'];
            writetable(subjectTab, tempfile);
            try
                exp = Experiment(tempfile, 'toydata', '');
            catch err
                delete(tempfile);
                rethrow(err);
            end
            delete(tempfile);

            exp.record.tabpath = fullfiletol(null_dir, 'toy_combined.csv');
            exp.name       = string('toy_null_vs_hypothesis');
            exp.subjectTab = subjectTab;

            fprintf('ToyDataLoader: loading %d null + %d hypothesis replicas...\n', ...
                    numel(null_paths), numel(hyp_paths));
            for i = 1:n
                exp.traces{i} = ToyDataLoader.loadTraces(all_paths{i}, all_groups{i});
                fprintf('  [%d/%d]  %s  (%s)\n', i, n, all_paths{i}, all_groups{i});
            end
            fprintf('Done.\n');
        end

    end

end

% ========================================================================= %
% Local helpers
% ========================================================================= %

function group = infer_group(s)
% INFER_GROUP  Derive a group label from the saved params when not provided.
    try
        if s.params.B == 0 && s.params.gamma == 0
            group = 'null';
        else
            group = 'hypothesis';
        end
    catch
        group = 'toy';
    end
end
