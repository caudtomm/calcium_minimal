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

        function manifolds = outputDataFiles(obj, outdir) % assumption: manifolds are 1 per stimulus per subject
            arguments
                obj GCMC_Analysis
                outdir char = 'manifold_data';
            end
            % prepare data for GCMC analysis
            % single fish
            % 
            v = obj.viewer;

            %% getting a list of the stimuli
            [~,labs] = v.dataFilter.filterData(v);
            stims = unique(labs{1});
            nStims = numel(stims);
            nSubjects = numel(labs);

            %% data extraction
            manifolds = cell(nSubjects,nStims);
            for i_stim = 1:nStims
                % only get this stimulus's responses
                v.dataFilter.stims_allowed = stims(i_stim);
                thispoints = v.dataFilter.filterData(v);

                % concatenate all repetitions (one stimulus = one manifold)
                thispoints = cellfun(@ActivityTraces.format, thispoints, 'UniformOutput',false);

                % eliminate any NaN
                idx_tokeep = cellfun(@(x) sum(isnan(x),2)==0, thispoints, 'UniformOutput',false);
                for i_sub = 1:nSubjects
                    thispoints{i_sub} = thispoints{i_sub}(idx_tokeep{i_sub},:);
                end

                % [t,N] -> [N,t]
                thispoints = cellfun(@transpose, thispoints, 'UniformOutput',false);

                % store to manifolds
                manifolds(:,i_stim) = thispoints;
            end

            %% save each subject's data to a python-compatible mat file
            if ~isfolder(outdir); mkdir(outdir); end
            for i_sub = 1:nSubjects
                data = manifolds(i_sub,:);
                save(fullfiletol(outdir,['manifolds_subj',num2str(i_sub),'.mat']), "data", "stims");
            end

        end

        function results = extractResults(obj,indir)
            arguments
                obj GCMC_Analysis
                indir char = pwd
            end

            v = obj.viewer;

            %% read and parse results

            % getting stimulus names ( # TODO : this should be retrieved from the output files,
            % which in turn should inherit it from their input file
            v.dataFilter.subjectIDs = {'TC_240104_TC0028_240101beh1b3_sxpDp_odorexp004_RPB3144501500AG'};
            v.dataFilter.stims_allowed = 'all stimuli';
            [~,labs] = v.dataFilter.filterData(v);
            stims = unique(labs{1});

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
                t = parseMetrics2table(allresults{i});
                
                % add the stimulus names
                t.manifold_name_1 = stims(t.manifold_idx_1 +1);
                t.manifold_name_2 = stims(t.manifold_idx_2 +1);

                results = [results; t];
            end

            %% store results
            outfname = fullfiletol(indir, 'GCMC_results_table.mat');
            save(outfname, 'results');
            disp(['Saved GCMC results table to ', outfname]);
        end

    end

end

%% helper functions

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