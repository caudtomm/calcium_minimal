function prettyJson = prettyjson(jsonData)
    % Create a version of JSON data with indentations and new lines

    % Initialize variables
    prettyJson = '';
    indent = 0;
    inString = false;
    
    % Loop through each character in the JSON data
    for i = 1:length(jsonData)
        char = jsonData(i);
        
        % Check if we are inside a string
        if char == '"'
            inString = ~inString;
        end
        
        % If not inside a string, format the JSON
        if ~inString
            switch char
                case '{'
                    indent = indent + 1;
                    prettyJson = [prettyJson, char, newline, repmat('    ', 1, indent)];
                case '}'
                    indent = indent - 1;
                    prettyJson = [prettyJson, newline, repmat('    ', 1, indent), char];
                case '['
                    indent = indent + 1;
                    prettyJson = [prettyJson, char, newline, repmat('    ', 1, indent)];
                case ']'
                    indent = indent - 1;
                    prettyJson = [prettyJson, newline, repmat('    ', 1, indent), char];
                case ','
                    prettyJson = [prettyJson, char, newline, repmat('    ', 1, indent)];
                otherwise
                    prettyJson = [prettyJson, char];
            end
        else
            prettyJson = [prettyJson, char];
        end
    end
end
