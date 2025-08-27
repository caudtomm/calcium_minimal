classdef ModeSelector
    properties
        v ExperimentViewer
        dataFilter DataFilter
        mode_name char = 'native_units' % Name of the mode
        values cell 
        coeffs cell % [subjects, 1] cell of double [units, ]
        data cell % [subjects, 1] cell of double [time, units, trials]
        labels cell
        params struct = struct() % Additional parameters (mode specific)
    end

    methods (Static)
        function values = calculateValues(coeffs, data)
            % Static method to calculate values based on coefficients and data
            values = cell(size(data));
            for i = 1:length(data)
                thisdata = data{i};
                [nTime, nUnits, nTrials] = size(thisdata);
                thisdata = ActivityTraces.format(thisdata); % [time*trials, units]
                thisvalues = thisdata * coeffs{i}; % # TODO might need fillmissing
                thisvalues = ActivityTraces.format(thisvalues, nTrials); % [time, modes, trials]
                values{i} = thisvalues;
            end
        end
    end
    
    methods
        function obj = ModeSelector(v, mode_name, dataFilter, varargin)
            arguments
                v ExperimentViewer
                mode_name char = []
                dataFilter = []
            end
            arguments (Repeating)
                varargin
            end

            % Parse single value inputs
            obj.v = v;
            if ~isempty(mode_name) && ischar(mode_name)
                obj.mode_name = mode_name; % input
            else
                obj.mode_name = 'native_units'; % default
            end
            if ~isempty(dataFilter) && isa(dataFilter,'DataFilter')
                obj.dataFilter = dataFilter; % input
            else
                obj.dataFilter = v.dataFilter; % default
            end
            if ~isempty(dataFilter); obj.dataFilter = dataFilter; end
            
            % Parse name-value pair arguments and store them in params
            obj.params = struct(varargin{:});

            % Specify input data
            [obj.data, obj.labels] = obj.dataFilter.filterData(v);
            
            % if pooling across subjects, concatenate data
            if isfield(obj.params, 'do_pool') && ~isempty(obj.params.do_pool)
                try
                    obj.data = {cat(2, obj.data{:})}; % concatenate along units dimension
                catch
                    warning('Could not pool data across subjects. Check dimensions. Proceeding without pooling.');
                end
            end

            % extract mode-specific values and coefficients
            obj = obj.extract;

        end
        
        function obj = extract(obj)
            % Extract mode-specific coefficients
            switch obj.mode_name
                case 'native_units'
                    obj.coeffs = cellfun(@(x) diag(ones(1,width(x))), obj.data, 'UniformOutput', false);
               case {'pca', 'nmf', 'ica', 'dpca'}
                    obj.coeffs = cell(size(obj.data));
                    for i = 1:length(obj.data)
                        thisdata = ActivityTraces.format(obj.data{i});
                        tempout = doModeDecomposition(...
                            thisdata, ...
                            'method', obj.mode_name, ...
                            'nfactors', obj.params.nfactors);
                        obj.coeffs{i} = tempout.coeffs;
                    end
                otherwise
                    error('Unknown mode name: %s', obj.mode_name);
            end

            obj.values = obj.calculateValues(obj.coeffs, obj.data);
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
    end
end