function outputpath = fullfiletol(thispath)
% tolerant fullfile, works with Unix and Windows

% Combine input arguments into a single path
combinedPath = fullfile(varargin{:});

% the external fullfile converts all Unix slashes to Windows slashes, if applicable.
outputpath = fullfile(strrep(combinedPath,'\','/'));

end
