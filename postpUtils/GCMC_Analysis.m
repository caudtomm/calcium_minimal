classdef GCMC_Analysis
    properties
        viewer ExperimentViewer
    end
    methods
        function obj = GCMC_Analysis(viewer)
            arguments
                viewer ExperimentViewer
            end
            obj.viewer = viewer;
        end

        function [pathlist, manifolds] = outputDataFiles_SlidingWindow(obj, mode, outdir, win_lim, win_size, win_step)
            arguments
                obj GCMC_Analysis
                mode char {mustBeMember(mode,{'odor','repetition','trial'})} = 'odor'
                outdir char = 'manifold_data'
                win_lim double = [0.5 20]
                win_size double = 4
                win_step double = 2
            end

            % Generate sliding windows
            win_start = win_lim(1):win_step:win_lim(2) - win_size;
            win_end = win_start + win_size;

            % Initialize output
            pathlist = strings(1, numel(win_start));
            %manifolds = 

            % Loop through each window
            for i = 1:numel(win_start)
                % Define the subdirectory for the current window
                subdir = fullfile(outdir, sprintf('%g-%gs', win_start(i), win_end(i)));

                % Call outputDataFiles() for the current window
                obj.viewer.dataFilter.interval = [win_start(i), win_end(i)];
                manifolds = obj.outputDataFiles(mode, subdir);

                % Store the subdirectory path
                pathlist(i) = subdir;
            end
        end

        function manifolds = outputDataFiles(obj, mode, outdir) % assumption: manifolds are 1 per stimulus per subject
            arguments
                obj GCMC_Analysis
                mode char {mustBeMember(mode,{'odor','repetition','trial'})} = 'odor'
                outdir char = 'manifold_data'
            end
            % prepare data for GCMC analysis
            % 

            v = obj.viewer;
            
            %% preparations (they assume all subjects have the same stimuli and trials)
            [~,labs] = v.dataFilter.filterData(v);
            stims = unique(labs{1});
            nTrials = numel(labs{1});
            reps = v.dataFilter.repetitions;
            if isempty(reps); reps = 1:sum(ismember(labs{1},stims(1))); end
            nSubjects = numel(labs);

            %% data extraction
            disp('Extracting manifolds for GCMC analysis...');
        
            % pick parameters based on mode
            switch mode
                case 'odor'
                    % all repetitions of this odor
                    % only get this stimulus's responses
                    filter_name = 'stims_allowed';
                    filter_vals = stims;
                case 'repetition'
                    % manifolds are repetitions across stimuli
                    filter_name = 'repetitions';
                    filter_vals = reps;
                case 'trial'
                    % manifolds are trials across stimuli (P = #trials * #stimuli)
                    filter_name = 'trial';
                    filter_vals = 1:nTrials;
            end

            manifolds = getManifolds(v, nSubjects, filter_name, filter_vals);

            %% save each subject's data to a python-compatible mat file
            disp(['Saving manifold data to ', outdir, '...']);
            if ~isfolder(outdir); mkdir(outdir); end
            [allSubjIDs, isinlist] = v.dataFilter.getSubjectIDs(v.subjectTab);
            subj_ordID = find(isinlist);
            for i_sub = 1:nSubjects
                data = manifolds(i_sub,:);
                subjID = allSubjIDs{i_sub};
                this_ordid = subj_ordID(i_sub);
                
                % for trial manifolds, store the whole stimulus id vector (one name per manifold!)
                if strcmp(mode,'trial')
                    stims = labs{i_sub};
                end
                
                save(fullfiletol(outdir,['manifolds_subj',num2str(this_ordid),'.mat']), "data", "stims", "subjID");
            end

            disp('Done.');

        end

        function [results, avg_results] = extractResults(obj,indir)
            arguments
                obj GCMC_Analysis
                indir char = pwd
            end

            v = obj.viewer;

            %% preparations
            
            % initialize outputs
            results = table;
            avg_results = table;

            % check that indir exists and contains mat files
            if ~isfolder(indir)
                warning(['Input directory does not exist: ', indir]);
                return;
            end
            files = dir(fullfiletol(indir,'*.mat'));
            if isempty(files)
                warning(['No .mat files found in input directory: ', indir]);
                return;
            end

            % prepare output file
            outfname = fullfiletol(indir, 'GCMC_results_table.mat');
            if isfile(outfname)
                disp(['Existing GCMC results table found: ', outfname,'. Overwriting...']);
                delete(outfname);
            end

            %% read and parse results

            % getting stimulus names ( # TODO : this should be retrieved from the output files,
            % which in turn should inherit it from their input file
            v.dataFilter.subjectIDs = {'TC_240104_TC0028_240101beh1b3_sxpDp_odorexp004_RPB3144501500AG'};
            v.dataFilter.stims_allowed = 'all stimuli';
            [~,labs] = v.dataFilter.filterData(v);
            oldstims = unique(labs{1}); % initialize

            files = dir(fullfiletol(indir,'*.mat'));

            % extract all results from individual files (pairwise manifold comparisons)
            allresults = cell(numel(files),1);
            for i = 1:numel(files)
                disp(['loading ',files(i).name])

                p = parseFileName(files(i).name);
                data = load(fullfiletol(indir,files(i).name));
                
                data.manifold_idx_1 = p.manifold_idx_1;
                data.manifold_idx_2 = p.manifold_idx_2;
                allresults{i} = data;
            end

            %% combine into a table
            results = table;
            for i = 1:numel(allresults)
                try
                    thisstims = allresults{i}.stims_i;
                    stims = thisstims(:);
                    oldstims = stims;
                catch
                    stims = oldstims;
                end
                t = parseMetrics2table(allresults{i});
                
                % add the stimulus names
                t.manifold_name_1 = stims(t.manifold_idx_1 +1);
                t.manifold_name_2 = stims(t.manifold_idx_2 +1);

                results = [results; t];
            end

            %% average over repetitions with the same parameters
            avg_results = genResultAvgTable(results, [], stims);

            %% store results
            save(outfname, 'results', 'avg_results');
            disp(['Saved GCMC results table to ', outfname]);
        end

        function [all_results] = extractResultsFromMultipleSubjects(obj, indir)
            arguments
                obj GCMC_Analysis
                indir char = pwd
            end

            files = dir(fullfiletol(indir,'manifolds_subj*.mat'));
            nSubjects = numel(files);

            all_results = struct();
            counter = 0;
            for i_sub = 1:nSubjects
                subj_name = ['manifolds_subj',num2str(i_sub)];
                subdir = fullfiletol(indir,subj_name);
                if exist(subdir,'dir')==0
                    disp(['Skipping subject ', num2str(i_sub), ': folder not found - ', subdir]);
                    continue;
                end
                counter = counter+1;
                disp(['Extracting results for subject ', num2str(i_sub), ' from ', subdir]);

                [results, avg_results] = obj.extractResults(subdir);
                all_results(counter).name = subj_name;
                all_results(counter).num_id = i_sub;
                all_results(counter).results = results;
                all_results(counter).avg_results = avg_results;

            end
        end

        function all_results = extractResults_SlidingWindow(obj,indir,subdir_name)
            arguments
                obj GCMC_Analysis
                indir char = pwd
                subdir_name char = ''
            end

            % Get a list of subfolders in the input directory
            subfolders = dir(indir);
            subfolders = subfolders([subfolders.isdir] & ~ismember({subfolders.name}, {'.', '..'}));

            % Initialize output structure
            all_results = struct();

            % Loop through each subfolder
            for i = 1:numel(subfolders)
                subfolder_path = fullfile(indir, subfolders(i).name, subdir_name);
                disp(['Processing subfolder: ', subfolder_path]);

                % Call extractResultsFromMultipleSubjects for the current subfolder
                results = obj.extractResultsFromMultipleSubjects(subfolder_path);

                % Store results in the output structure
                all_results(i).data = results;
                all_results(i).window = subfolders(i).name;
            end

            % Reorder according to naming, if possible
            all_results = sortByWindow(all_results);
        end

        function group_data = clusterByGroup(obj, all_results)
            arguments
                obj GCMC_Analysis
                all_results struct
            end
            
            metrics = all_results(1).avg_results.Properties.VariableNames(5:end-2); % exclude grouping variables and stimulus names
            subjectTab = obj.viewer.subjectTab;
            nResults = numel(all_results);
            
            % select subjectTab rows (subjects) according to the results'
            % ordinal ID
            res_ids = [];
            for i = 1:nResults
                this_id = all_results(i).num_id;
                res_ids = [res_ids; this_id];
            end
            subjectTab = subjectTab(res_ids,:);

            % separate data by subject group
            groups = unique(subjectTab.group);
            nGroups = numel(groups);
            group_data = struct();
            for i = 1:nGroups
                idx = ismember(subjectTab.group,groups{i});
                this_res = all_results(idx);
                group_data(i).group_name = groups{i};
                group_data(i).data = table();

                for i_subj = 1:numel(this_res)
                    this_subj_id = this_res(i_subj).num_id;
                    this_data = this_res(i_subj).avg_results;
                    if isempty(this_data); continue; end % escape
                    this_data.subj_id = repmat(this_subj_id, height(this_data),1);
                    group_data(i).data = [group_data(i).data; this_data];
                end
            end
        end

    end

    methods (Static)
        function new_group_data = catTrainedGroups(group_data)
            % concat trained groups (brittle, but I don't really need more..)
            new_group_data = struct;
            new_group_data(1).group_name = group_data(1).group_name;
            new_group_data(1).data = group_data(1).data;
            new_group_data(2).group_name = 'trained';
            new_group_data(2).data = [group_data(2).data;group_data(3).data;group_data(4).data];
            new_group_data(3).data = group_data(5).data;
            new_group_data(3).group_name = group_data(5).group_name;
        end
    end
end

%% helper functions


function avg_results = genResultAvgTable(results, idx, stims)
    if nargin < 2 || isempty(idx)
        avg_results = varfun(@(x) mean(x,"omitmissing"), results, ...
            'InputVariables', results.Properties.VariableNames(6:end-2), ...
            'GroupingVariables', {'manifold_idx_1','manifold_idx_2','shuffle'});
    else
        avg_results = varfun(@(x) mean(x(idx,:),"omitmissing"), results, ...
            'InputVariables', results.Properties.VariableNames(6:end-2), ...
            'GroupingVariables', {'manifold_idx_1','manifold_idx_2','shuffle'});
    end

    % add the stimulus names
    if ~exist("stims","var") || isempty(stims); return; end
    avg_results.manifold_name_1 = stims(avg_results.manifold_idx_1 +1);
    avg_results.manifold_name_2 = stims(avg_results.manifold_idx_2 +1);
end

function T = parseMetrics2table(s)
    % s : results structure

    M = s.metrics;

    % front-append index names tp column names
    cols = string(M.columns);     % 1×K string array (column labels)
    idx_names = {'manifold_idx_1','manifold_idx_2','repetition','shuffle','seed'};
    cols = [idx_names, cols];

    % front-append index values to data matrix
    idx  = M.index;               % Nx1 (char cell) - format {'(rep,shuffle,seed)'}
    [rep,shuf,seed] = cellfun(@(x) parseIndexVal(x), idx, 'UniformOutput', false);
    rep = cell2mat(rep); shuf = cell2mat(shuf); seed = cell2mat(seed);
    N = numel(rep);
    mid1 = repelem(s.manifold_idx_1,N,1); mid2 = repelem(s.manifold_idx_2,N,1);
    X    = M.values;              % NxK double
    X    = [mid1,mid2,rep,shuf,seed, X]; % Nx(K+7) double
    
    T = array2table(X, 'VariableNames', cellstr(cols));
    
end

function [firstInt, boolVal, lastInt] = parseIndexVal(s)
    % s : char vector or string
    
    % Extract integers
    nums = regexp(s, '\d+', 'match');
    firstInt = str2double(nums{1});
    lastInt  = str2double(nums{end});
    
    % Extract boolean
    boolStr = regexp(s, '(True|False)', 'match', 'ignorecase');
    boolVal = strcmpi(boolStr{1}, 'True');
end

function out = parseFileName(fname)
    % expected format: p_<manifold_idx_1>_<manifold_idx_2>.mat
    tokens = regexp(fname, '^p_(\d+)_(\d+)\.mat$', 'tokens');
    if isempty(tokens)
        error('Filename does not match expected format: %s', fname);
    end

    % extract indices
    tokens = tokens{1};
    out.manifold_idx_1 = str2double(tokens{1});
    out.manifold_idx_2 = str2double(tokens{2});
end

function manifolds = getManifolds(v, nSubjects, filter_name, filter_vals) % this could be a method of ExperimentViewer
    % get manifolds grouped by repetition number (irrespective of stimulus identity)
    nVals = numel(filter_vals);

    % initialize output
    manifolds = cell(nSubjects,nVals);

    switch filter_name
    case 'trial'
        [~, thispoints] = ModeSelector(v).extract;

        % loop over each individual trial number and subject
        for i_subj = 1:nSubjects
            subj_points = thispoints{i_subj};
            if isempty(subj_points); continue; end
            for i_val = 1:nVals
                % extract this trial's points for this subject
                trial_points = {subj_points(:,:,i_val)}; % {[t,N] double}

                % store to manifolds
                manifolds(i_subj,i_val) = formatManifold(trial_points);
            end
        end
        
    otherwise
        for i_val = 1:nVals
            v.dataFilter.(filter_name) = filter_vals(i_val);
            
            [~, thispoints] = ModeSelector(v).extract;

            % store to manifolds
            manifolds(:,i_val) = formatManifold(thispoints);
        end
    end
end

function manifold = formatManifold(points)
    % points : cell array of size nSubjects x 1, each cell is [t,N] double
    nSubjects = numel(points);

    % initialize output
    manifold = {};

    % check for empty input
    if isempty(points); return; end

    % concatenate all observations (one filter value = one manifold)
    try
        points = cellfun(@ActivityTraces.format, points, 'UniformOutput',false);
    catch
    end

    % eliminate any NaN
    idx_tokeep = cellfun(@(x) sum(isnan(x),2)==0, points, 'UniformOutput',false);
    for i_sub = 1:nSubjects
        points{i_sub} = points{i_sub}(idx_tokeep{i_sub},:);
    end

    % [t,N] -> [N,t]
    points = cellfun(@transpose, points, 'UniformOutput',false);

    % return
    manifold = points;
end

function [results_out, idx] = sortByWindow(results_in)

% initialize output
results_out = results_in;
if isempty(results_in); return; end

% Extract the numerical values from the .name field
names = {results_in.window};
names = strrep(names, '_', '.'); % for decimals
pattern = '([-]?\d+\.?\d*)-([-]?\d+\.?\d*)\w*';
numeric_values = nan(numel(names), 1);

for i = 1:numel(names)
    tokens = regexp(names{i}, pattern, 'tokens');
    if ~isempty(tokens)
        numeric_values(i) = str2double(tokens{1}{1}); % Extract the first %g
    end
end

% Sort based on the extracted numerical values
[~, idx] = sort(numeric_values, 'ascend');
results_out = results_in(idx);
end