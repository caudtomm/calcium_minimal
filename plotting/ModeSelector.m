classdef ModeSelector
    properties
        v ExperimentViewer
        mode_name char = 'native_units' % Name of the mode
        mode_method char = 'mode_values' % {'mode_values', 'isolate', 'subtract'}
        mode_OI = 'all' % modes of interest: {'all', 'stimulus', 'non-stimulus', 'novelty'}
        mode_file char = ''; % if empty, extracts by default. Else, it looks for coefficients in the file specified.
        values cell 
        coeffs cell % [subjects, 1] cell of double [units, ]
        data cell % [subjects, 1] cell of double [time, units, trials]
        labels cell
        params struct = struct() % Additional parameters (mode specific)

        fullout
    end

    methods (Static)
        function values = calculateValues(coeffs, data)
            % Static method to calculate values based on coefficients and data
            values = cell(size(data));
            for i = 1:length(data)
                thisdata = data{i};
                [nTime, nUnits, nTrials] = size(thisdata);
                thisdata = ActivityTraces.format(thisdata); % [time*trials, units]
                if isempty(thisdata); continue; end
                thisvalues = thisdata * coeffs{i};
                thisvalues = ActivityTraces.format(thisvalues, nTrials); % [time, modes, trials]
                values{i} = thisvalues;
            end
        end
    end
    
    methods
        function obj = ModeSelector(v, varargin)
            arguments
                v ExperimentViewer
            end
            arguments (Repeating)
                varargin
            end

            % Parse single value inputs
            obj.v = v;

            % Parse DataFilter-derived mode selection options
            dft = v.dataFilter;
            obj.mode_name = dft.mode_name;
            obj.mode_method = dft.mode_method;
            obj.mode_OI = dft.mode_OI;
            obj.mode_file = dft.mode_file;            
            obj.params = dft.mode_params;

            % Specify input data
            [obj.data, obj.labels] = dft.filterData(v);
            
            % if pooling across subjects, concatenate data
            if isfield(obj.params, 'do_pool') && ~isempty(obj.params.do_pool)
                try
                    obj.data = {cat(2, obj.data{:})}; % concatenate along units dimension
                catch
                    warning('Could not pool data across subjects. Check dimensions. Proceeding without pooling.');
                end
            end

        end

        % extract modes according to specified method
        function [obj, values, labels, coeffs] = extract(obj, force_import, force_recompute)
            arguments
                obj ModeSelector
                force_import logical = false
                force_recompute logical = false
            end

            if force_recompute; obj = obj.wipeResults; end

            % initialize output
            labels = obj.labels; % this is it for the labels
            values = obj.data; % compute values from scratch to avoid recursion
            coeffs = obj.coeffs;

            % options
            method = obj.mode_method;
            fname = obj.mode_file;
            nsubjects = numel(obj.data);

            [obj, coeffs] = obj.importWeightsFile(force_import);

            % compute if needed
            if isempty(coeffs) || force_recompute
                obj = obj.compute;
                coeffs = obj.coeffs;
            end
            
            values = obj.calculateValues(coeffs, obj.data);
            obj.values = values;
            if all(cellfun(@isempty,values)); return; end
            
            switch method
                case 'mode_values'
                    % only retain values of interest
                    for i = 1:nsubjects
                        mids = obj.parseModeOI(i);
                        values{i} = values{i}(:, mids, :);
                    end
                case 'isolate'
                    values = obj.isolate;
                case 'subtract'
                    values = obj.isolate;
                    % elementwise subtraction from obj.data
                    for i = 1:nsubjects
                        values{i} = obj.data{i} - values{i};
                    end
                otherwise
                    error('Unknown mode_method: %s', method);
            end

            obj.values = values;
        end
        
        % compute mode-specific coefficients
        function obj = compute(obj)
            disp(['Computing coefficients. Mode: ',obj.mode_name])
            % Extract mode-specific coefficients
            switch obj.mode_name
                case 'native_units'
                    obj.coeffs = cellfun(@(x) diag(ones(1,width(x))), obj.data, 'UniformOutput', false);
                case {'pca', 'nmf', 'ica'}
                    obj.coeffs = cell(size(obj.data));
                    for i = 1:length(obj.data)
                        thisdata = ActivityTraces.format(obj.data{i});
                        if ~isfield(obj.params,'nfactors')
                            % use defaults
                            tempout = doModeDecomposition(...
                                thisdata, ...
                                'method', obj.mode_name);
                        else
                            tempout = doModeDecomposition(...
                                thisdata, ...
                                'method', obj.mode_name, ...
                                'nfactors', obj.params.nfactors);
                        end
                        obj.coeffs{i} = tempout.coeffs;
                    end
                case 'dpca'
                    for i = 1:length(obj.data)
                        nReps = obj.params.nReps;
                        framerate = obj.params.framerate;

                        % prepare data for dpca
                        thisdata = obj.data{i}; % [time, units, trials]
                        [nTime, nUnits, nTrials] = size(thisdata);
                        thisdata = ActivityTraces.format(thisdata);
                        thisdata = fillmissing(thisdata,"previous");
                        thisdata = fillmissing(thisdata,"constant",0);
                        thisdata = ActivityTraces.format(thisdata,nTrials);
                        thislabels = obj.labels{i};
                        [thislabels, idx] = sort(thislabels); % sort labels
                        thisdata = thisdata(:,:,idx); % sort data accordingly
                        uniqueLabels = unique(thislabels, 'rows');
                        nLabels = size(uniqueLabels, 1);
                        tmp = nan(nTime, nUnits, nLabels, nReps);
                        for l = 1:nLabels
                            labelMask = ismember(thislabels, uniqueLabels(l,:), 'rows');
                            tmp(:,:,l,:) = thisdata(:,:,labelMask);
                        end
                        thisdata = tmp; % [time, units, labels, repetitions]

                        % prepare inputs for dpca
                        trialNum = nReps * ones(nUnits,nLabels,1); % [N, S, D]
                        firingRates = permute(thisdata,[2,3,5,1,4]); % [N, S, D, T, E]
                        firingRatesAverage = mean(firingRates,5,'omitmissing'); % [N, S, D, T]
                        t = (1:nTime)/framerate;

                        % run dpca
                        out = dpca_fromdemo(trialNum, firingRates, firingRatesAverage, t, 0);
                        
                        % store output
                        obj.coeffs{i} = out.W;
                        obj.fullout{i} = out;
                    end

                otherwise
                    error('Unknown mode name: %s', obj.mode_name);
            end

        end
        
        % recontruct unit activity based only on a subset of modes
        function values = isolate(obj)
            values = obj.values;
            coeffs = obj.coeffs;
            nsubjects = numel(obj.data);

            for i_sj = 1:nsubjects
                thiscoeffs = coeffs{i_sj};
                thisvalues = values{i_sj};
                [T, nmodes, nTrials] = size(thisvalues);
                thisvalues = ActivityTraces.format(thisvalues); % [T * nTrials x nmodes]

                % specify ids for modes of interest
                mids = obj.parseModeOI(i_sj);

                % zero-out values for any unselected modes
                zeromat = zeros(T*nTrials, sum(~mids));
                thisvalues(:, ~mids) = zeromat;
            
                % get inverse weight matrix
                W = pinv(thiscoeffs); % [modes x units]

                % reconstruct unit activity
                thisvalues = thisvalues * W; % [T * nTrials x N]

                % refold 
                values{i_sj} = ActivityTraces.format(thisvalues, nTrials); % [T x N x nTrials]
            end
        end

        function mids = parseModeOI(obj, subjectNum)
            mode_OI = obj.mode_OI;
            thiscoeffs = obj.coeffs{subjectNum};

            if ischar(mode_OI)
            mids = obj.getModeIDs(subjectNum); % logical
            elseif isnumeric(mode_OI)
            mids = false(size(thiscoeffs, 2), 1);
            mids(mode_OI) = true;
            elseif islogical(mode_OI)
            mids = mode_OI;
            if length(mids) ~= size(thiscoeffs, 2)
                error('Length of logical mode_OI does not match number of modes.');
            end
            elseif iscell(mode_OI)
            if subjectNum > numel(mode_OI)
                error('Subject number exceeds the number of elements in mode_OI cell array.');
            end
            mids = obj.parseModeOI(subjectNum); % Recursive call on the content of mode_OI{subjectNum}
            else
            error('mode_OI must be char, numeric, logical, or cell.');
            end
        end

        % turn mode_OI string into logical mode ids
        function idx = getModeIDs(obj, subjectNum)
            [N,nmodes] = size(obj.coeffs{subjectNum});
            moistr = obj.mode_OI;
            if numel(obj.fullout)>=subjectNum && isfield(obj.fullout{subjectNum}, 'whichMarg')
                whichMarg = obj.fullout{subjectNum}.whichMarg;
                stimulus_marginalizations = ismember(whichMarg, [1,2,4]);
            elseif ~strcmp(moistr, 'all')
                warning('Field "whichMarg" does not exist in fullout{%d}. Defaulting mode_OI to "all".', subjectNum);
                moistr = 'all';
            end

            idx = true(nmodes,1);

            switch moistr
                case 'all'
                    % do nothing
                case 'all_stimulus'
                    idx = stimulus_marginalizations;
                case 'stimulus'
                    idx = stimulus_marginalizations;
                    idx(find(idx,1)) = false; % eliminate the highest variant one
                case 'non-stimulus'
                    idx = ~stimulus_marginalizations;
                case 'novelty'
                    idx = stimulus_marginalizations;
                    idx(find(idx,1)) = false; % eliminate the highest variant one
                    idx = ~idx; % keep only the highest variant stimulus mode ("novelty")
                otherwise
                    warning('Unknown mode_OI: %s. Returning all modes.', moistr);
            end
        end

        function [values, coeffs] = binarize(obj, threshold)
            arguments
                obj ModeSelector
                threshold double = 0
            end

            % initialize output
            values = obj.data;
            coeffs = obj.coeffs;

            % binarize coeffs
            coeffs = coeffs >= threshold;

            % recalculate values
            values = obj.calculateValues(coeffs);

            
        end

        function [obj, coeffs] = importWeightsFile(obj, force_import)
            % Initialize output
            coeffs = obj.coeffs;
            fname = obj.mode_file;

            % If no file specified, skip
            if isempty(fname)
                disp('No weights file specified. Skipping import.');
                return;
            end

            % Check if file exists and load weights
            if isempty(coeffs) || force_import
            if exist(fname, "file")
                try
                disp(['Loading weights from file: ', fname]);
                fileIn = load(fname);
                
                % Check if the required field exists in the loaded file
                if isfield(fileIn, 'dpca')
                    obj.fullout = fileIn.dpca;
                    coeffs = cellfun(@(x) x.W, fileIn.dpca, 'UniformOutput', false);
                    obj.coeffs = coeffs;
                    disp('Weights successfully loaded and assigned.');
                else
                    warning('The file does not contain the expected "dpca" field. Skipping weight import.');
                end
                catch ME
                warning(['An error occurred while loading weights from file: ', fname]);
                disp(['Error message: ', ME.message]);
                end
            else
                warning(['Specified weights file does not exist: ', fname]);
            end
            else
            disp('Coefficients are already initialized. Skipping file import.');
            end
        end
        
        function obj = wipeResults(obj)
            obj.values = {};
            obj.coeffs = {};
            obj.fullout = {};
        end

        % plotting by external function hf = genFigures(obj)
    end

end