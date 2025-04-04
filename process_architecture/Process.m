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
            hashlen = 20;
            hash = char(randi([97 122],1,hashlen));
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

        function process = getLight(obj)
            % empties out all properties considered memory-heavy
            arguments
                obj
            end
            process = obj;
            process = process.clearProperties(process.heavyProperties);
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