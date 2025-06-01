classdef (Abstract) Process
    properties
        init ProcessInit
    end
    properties (SetAccess = protected)
        hash char
    end
    properties (SetAccess = private)
        heavyProperties = { 'data_raw';
                            'data_processed'}
    end

    
    methods (Abstract, Static)
        run
    end

    methods (Abstract, Access=protected)
        disp_runheader
    end

    methods (Static)
        function hash = setHash()
            rng('shuffle');  % Seed the random number generator based on current time
            hashlen = 20;
            hash = char(randi([97 122],1,hashlen));

            % % alternative generation, more robust for unique hashes, but
            % % can contain '-' and a couple other special characters
            % hashlen = 20;
            % uuid = char(java.util.UUID.randomUUID.toString());
            % hash = uuid(1:hashlen);  % Trim to 20 characters if needed
        end
    end

    methods
        function obj = Process() % Constructor
            obj.init = ProcessInit;
            obj.init.processclass = class(obj);
            obj = obj.resetHash;
        end

        function obj = resetHash(obj)
            obj.hash = obj.setHash;
            obj.init.hash = obj.hash;
        end

        function process = save(obj, filename)
            arguments
                obj Process
                filename char = [obj.hash,'.mat']
            end
            if isfolder(filename)
                filename = fullfiletol(filename, [obj.hash,'.mat']);
            end
            process = obj;
            s.process = process;
            robust_io('save',filename,s,'-mat','-v7.3')
        end

        function process = getLight(obj, recursive)
            % Recursively empties out all properties considered memory-heavy
            arguments
                obj
                recursive logical = false
            end
            process = obj;
            if recursive
                process = recursiveClear(process);
            else
                process = process.clearProperties(process.heavyProperties);
            end
        end
        
        function obj = recursiveClear(obj)
            % Clear current object's heavy properties
            obj = obj.clearProperties(obj.heavyProperties);
        
            % Loop over all properties of the object
            mc = metaclass(obj);
            props = mc.PropertyList;
            for i = 1:numel(props)
                field = props(i).Name;
                value = obj.(field);

                if strcmp(props(i).SetAccess,'private'); continue; end
        
                % If the property is the same class and not empty, apply recursively
                if isa(value, class(obj)) && ~isempty(value)
                    obj.(field) = recursiveClear(value);
                elseif iscell(value)
                    for j = 1:numel(value)
                        if isa(value{j}, class(obj)) && ~isempty(value{j})
                            value{j} = recursiveClear(value{j});
                        end
                    end
                    obj.(field) = value;
                elseif isstruct(value)
                    for j = 1:numel(value)
                        value(j) = structRecursiveClear(value(j), class(obj));
                    end
                    obj.(field) = value;
                end
            end
        end
        
        function s = structRecursiveClear(s, className)
            fns = fieldnames(s);
            for k = 1:numel(fns)
                fn = fns{k};
                val = s.(fn);
                if isa(val, className) && ~isempty(val)
                    s.(fn) = recursiveClear(val);
                elseif iscell(val)
                    for j = 1:numel(val)
                        if isa(val{j}, className) && ~isempty(val{j})
                            val{j} = recursiveClear(val{j});
                        end
                    end
                    s.(fn) = val;
                elseif isstruct(val)
                    for j = 1:numel(val)
                        val(j) = structRecursiveClear(val(j), className);
                    end
                    s.(fn) = val;
                end
            end
        end


        function obj = clearProperties(obj, propNames)
            % clearProperties Clears the contents of specified properties of an object
            % while maintaining their original types.
            % 
            % Inputs:
            %   obj - The object whose properties are to be modified.
            %   propNames - A cell array of strings specifying the property names to clear.
            
            % Get the meta-object to access property information
            metaObj = metaclass(obj);
            % Get list of fields in obj
            objfields = {metaObj.PropertyList.Name};
            % Iterate over the list of property names
            for k = 1:length(propNames)
                propName = propNames{k};
                % Check if the property exists in the object
                if any(strcmp(objfields, propName))
                    % Depending on the property type, clear its content accordingly
                    switch class(obj.(propName))
                        case 'double'
                            obj.(propName) = double.empty;
                        case 'char'
                            obj.(propName) = '';
                        case 'cell'
                            obj.(propName) = {};
                        case 'struct'
                            obj.(propName) = struct();
                        case 'Movie'
                            obj.(propName).stack = [];
                            obj.(propName).scanimage_meta = {};
                        % Add cases for other types as needed
                        otherwise
                            % Handle generic case or unknown types
                            warning('Type of %s is not explicitly handled. Setting to empty.', propName);
                            obj.(propName) = [];
                    end
                else
                    % Property does not exist; issue a warning
                    warning('Property %s does not exist in object of class %s.', propName, class(obj));
                end
            end

            % recursion
            for i_field = 1:numel(objfields)
                fieldName = objfields{i_field};
                if isSubclassOf(obj.(fieldName),'Process')
                    try; obj.(fieldName) = getLight(obj.(propName)); end
                end
            end

        end

    end
end