function propStruct = objectToPropsStruct(obj, excludeProps)
    % objectToPropsStruct Converts object properties to a structure, excluding specified properties.
    % 
    % Inputs:
    %   obj - The object whose properties are to be converted.
    %   excludeProps - A cell array of strings listing the properties to exclude.

    % Get meta information about the object's class
    metaObj = metaclass(obj);
    allProps = metaObj.PropertyList;

    % Initialize an empty structure
    propStruct = struct();

    % Loop through all properties
    for k = 1:length(allProps)
        propName = allProps(k).Name;
        
        % Check if the current property should be excluded
        if ~ismember(propName, excludeProps)
            % Only add the property if it's not in the exclude list
            propValue = obj.(propName);
            propStruct.(propName) = propValue;
        end
    end
end
