function varargout = robust_io(mode, filepath, varargin)
%ROBUST_IO Load or save with retry logic
%
% Usage:
%   data = robust_io('load', filepath)
%   robust_io('save', filepath, dataStruct)
%
% Automatically retries on failure (e.g. flaky network I/O)
%
% For saving, pass a struct containing variables to save:
%   s.myVar = myVar;
%   s.other = otherVar;
%   robust_io('save', filepath, s);

    maxRetries = 5;
    delay = 30;  % seconds
    attempt = 1;
    success = false;

    while ~success && attempt <= maxRetries
        try
            switch lower(mode)
                case 'load'
                    data = load(filepath, varargin{:});
                    varargout{1} = data;
                case 'save'
                    if isempty(varargin)
                        error('No variables provided to save.');
                    end
                    varsToSave = varargin{1};
                    if ~isstruct(varsToSave)
                        error('For saving, provide a struct of variables to save.');
                    end
                    save(filepath, '-struct', 'varsToSave', varargin{2:end});
                otherwise
                    error('Unsupported mode: %s', mode);
            end
            success = true;
        catch ME
            fprintf('I/O attempt %d failed: %s\n', attempt, ME.message);
            if attempt == maxRetries
                rethrow(ME);
            else
                pause(delay);
                attempt = attempt + 1;
            end
        end
    end
end
