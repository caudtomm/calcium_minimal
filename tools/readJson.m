function outstruct = readJson(FileIn)
    fname = getFileNameSpecs(FileIn).fname;
    if ~exist(FileIn,'file'); warn_notfound(fname); return; end
    jsonStr = fileread(FileIn);
    try
        outstruct = jsondecode(jsonStr);
        showcontent(outstruct)
    catch
        warning('couldn''t read: %s', fname)
    end
end