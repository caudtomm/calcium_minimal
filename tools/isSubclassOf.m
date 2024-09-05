function isSubclass = isSubclassOf(obj, className)
    % Check if the class of obj is a subclass of the specified className.
    % Inputs:
    %   obj - the object to check.
    %   className - the name of the potential superclass.
    % Output:
    %   isSubclass - true if obj is a subclass of className, otherwise false.
    
    superClassNames = superclasses(obj);
    isSubclass = any(strcmp(superClassNames, className));
end
