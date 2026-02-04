classdef GCMCResultsFilter
    % GCMCResultsFilter - Flexible filtering for GCMC results tables
    %
    % Filters can be value-based (match specific values) or relational
    % (compare column pairs like "same name" or "different repetition").
    %
    % Example usage:
    %   % Basic filter
    %   f = GCMCResultsFilter('shuffle', false);
    %
    %   % Filter by specific stimuli with same-name pairs only
    %   f = GCMCResultsFilter('manifold_name_1', {'Arg', 'Ala'}, ...
    %                         'manifold_name_2', {'Arg', 'Ala'});
    %   f.pairFilters.manifold_name = 'same';  % Arg-Arg, Ala-Ala only
    %
    %   % Different-name pairs only
    %   f = GCMCResultsFilter();
    %   f.pairFilters.manifold_name = 'different';
    %
    %   % Apply filter
    %   filtered_table = f.filterTable(results_table);

    properties
        % Value-based filters (empty = no filter / include all)
        shuffle logical = []              % true/false to filter by shuffle status
        manifold_name_1 cell = {}         % filter by manifold_name_1 values
        manifold_name_2 cell = {}         % filter by manifold_name_2 values
        subj_ids double = []              % filter by subj_id values
        manifold_rep1 double = []         % filter by manifold_rep1 values
        manifold_rep2 double = []         % filter by manifold_rep2 values

        % Dynamic/relational filters (compare column pairs)
        % Format: struct where fieldname = column base name (e.g., 'manifold_name')
        %         value = 'same' or 'different'
        % Example: pairFilters.manifold_name = 'same' keeps rows where name_1 == name_2
        pairFilters struct = struct()

        % Generic column filters (for extensibility)
        % Format: struct where fieldname = exact column name, value = allowed values
        customFilters struct = struct()
    end

    methods
        function obj = GCMCResultsFilter(varargin)
            % Constructor: accepts struct or name-value pairs
            %   GCMCResultsFilter()
            %   GCMCResultsFilter(struct)
            %   GCMCResultsFilter('propName', value, ...)

            if nargin == 1 && isstruct(varargin{1})
                % Initialize from struct
                s = varargin{1};
                fn = fieldnames(s);
                for i = 1:numel(fn)
                    if isprop(obj, fn{i})
                        obj.(fn{i}) = s.(fn{i});
                    end
                end
            elseif mod(nargin, 2) == 0
                % Name-value pair input
                for i = 1:2:nargin
                    name = varargin{i};
                    value = varargin{i+1};
                    if isprop(obj, name)
                        obj.(name) = value;
                    else
                        error('Invalid property name: %s', name);
                    end
                end
            elseif nargin > 0
                error('Unsupported GCMCResultsFilter constructor usage.');
            end
        end

        function filtered = filterTable(obj, T)
            % Apply all filters and return filtered table
            arguments
                obj GCMCResultsFilter
                T table
            end

            idx = obj.getFilterMask(T);
            filtered = T(idx, :);
        end

        function idx = getFilterMask(obj, T)
            % Returns logical mask for rows to keep
            arguments
                obj GCMCResultsFilter
                T table
            end

            idx = true(height(T), 1);
            colNames = T.Properties.VariableNames;

            % Value-based filters
            if ~isempty(obj.shuffle)
                idx = idx & (T.shuffle == obj.shuffle);
            end

            if ~isempty(obj.manifold_name_1) && ismember('manifold_name_1', colNames)
                idx = idx & ismember(T.manifold_name_1, obj.manifold_name_1);
            end

            if ~isempty(obj.manifold_name_2) && ismember('manifold_name_2', colNames)
                idx = idx & ismember(T.manifold_name_2, obj.manifold_name_2);
            end

            if ~isempty(obj.subj_ids) && ismember('subj_id', colNames)
                idx = idx & ismember(T.subj_id, obj.subj_ids);
            end

            if ~isempty(obj.manifold_rep1) && ismember('manifold_rep1', colNames)
                idx = idx & ismember(T.manifold_rep1, obj.manifold_rep1);
            end

            if ~isempty(obj.manifold_rep2) && ismember('manifold_rep2', colNames)
                idx = idx & ismember(T.manifold_rep2, obj.manifold_rep2);
            end

            % Dynamic pair filters
            idx = idx & obj.applyPairFilters(T);

            % Custom filters
            cols = fieldnames(obj.customFilters);
            for i = 1:numel(cols)
                col = cols{i};
                if ismember(col, colNames)
                    idx = idx & ismember(T.(col), obj.customFilters.(col));
                end
            end
        end
    end

    methods (Access = private)
        function idx = applyPairFilters(obj, T)
            % Apply relational filters that compare column pairs
            idx = true(height(T), 1);
            cols = fieldnames(obj.pairFilters);
            colNames = T.Properties.VariableNames;

            for i = 1:numel(cols)
                base = cols{i};  % e.g., 'manifold_name' or 'manifold_rep'
                mode = obj.pairFilters.(base);  % 'same' or 'different'

                % Derive column names (try both naming conventions)
                col1 = [base '_1'];
                col2 = [base '_2'];
                if ~ismember(col1, colNames)
                    col1 = [base '1'];  % e.g., manifold_rep1
                    col2 = [base '2'];
                end

                if ismember(col1, colNames) && ismember(col2, colNames)
                    switch mode
                        case 'same'
                            idx = idx & compareColumns(T.(col1), T.(col2), @eq);
                        case 'different'
                            idx = idx & compareColumns(T.(col1), T.(col2), @ne);
                        otherwise
                            warning('Unknown pairFilter mode: %s. Ignoring.', mode);
                    end
                end
            end
        end
    end
end

%% Helper functions

function match = compareColumns(col1, col2, op)
    % Compare two columns element-wise, handling cell arrays and numeric
    if iscell(col1)
        if isequal(op, @eq)
            match = strcmp(col1, col2);
        else
            match = ~strcmp(col1, col2);
        end
    else
        match = op(col1, col2);
    end
end
