function [fname, fext, fpath] = getFName(str)
% gets part of string after the last '\' and before the last '.'.

% initialization
fname = deal(str);
fpath = '';

% fname + fext
fext = flip(extractBefore(flip(str),'.'));
if isempty(fext); warning('Extension missing from Filename'); end
fname = extractBefore(str, strcat('.',fext));

% extract file path if present
if ~isempty(extractBefore(fname, '\'))
    fpath = flip(extractAfter(flip(fname),'\'));
    fname = extractAfter(fname, strcat(fpath,'\'));
elseif ~isempty(extractBefore(fname, '/'))
    fpath = flip(extractAfter(flip(fname),'/'));
    fname = extractAfter(fname, strcat(fpath,'/'));
else
    warning('File path not found: left blank.')
    return
end

end