function fullfiletol(thispath)
% tolerant fullfile, works with Unix and Windows

% the internal fullfile acts normally to bind char vectors.
% the external fullfile converts all Unix slashes to Windows slashes, if applicable.
fullfile(strrep(fullfile(thispath),'\','/'))

end
