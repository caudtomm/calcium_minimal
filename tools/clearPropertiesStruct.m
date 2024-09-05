function s = clearPropertiesStruct(s, fieldNames)
    % clearPropertiesStruct Clears the contents of specified fields in a struct
    % while maintaining their original types.
    % 
    % Inputs:
    %   s - The struct whose fields are to be cleared.
    %   fieldNames - A cell array of strings specifying the field names to clear.
    
    % Get list of fields in the struct
    structFields = fieldnames(s);
    
    % Iterate over the list of field names
    for k = 1:length(fieldNames)
        fieldName = fieldNames{k};
        
        % Check if the field exists in the struct
        if any(strcmp(structFields, fieldName))
            % Depending on the field type, clear its content accordingly
            switch class(s.(fieldName))
                case 'double'
                    s.(fieldName) = double.empty;
                case 'char'
                    s.(fieldName) = '';
                case 'cell'
                    s.(fieldName) = {};
                case 'struct'
                    s.(fieldName) = struct();
                case 'Movie'
                    s.(fieldName).stack = [];
                otherwise
                    % Handle generic case or unknown types
                    warning('Type of %s is not explicitly handled. Setting to empty.', fieldName);
                    s.(fieldName) = [];
            end
        else
            % Field does not exist; issue a warning
            warning('Field %s does not exist in the struct.', fieldName);
        end
    end

    % Recursion: check if any field is a struct and clear its fields recursively
    for iField = 1:numel(structFields)
        currentField = structFields{iField};
        if isstruct(s.(currentField))
            s.(currentField) = clearPropertiesStruct(s.(currentField), fieldNames);
        end
    end
end