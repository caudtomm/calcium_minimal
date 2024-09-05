function daughterfile = find_daughter_file(basename,obligate_ext)
    % change cd
    thisdir = pwd;
    srcdir = getFileNameSpecs(basename).orig_fpath;
    if ~isempty(srcdir); cd(srcdir); end
    basename = getFileNameSpecs(basename).fname;

    % search is carried out in the current directory
    if ~isempty(obligate_ext)
        files = dir(['*.',obligate_ext]);
    else 
        files = dir();
    end
    files = {files(:).name}';

    % remove extension (if applicable)
    if contains(basename,'.'); basename = extractBefore(basename,'.'); end

    % search daughters
    idx = contains(files, basename);
    daughterfile = fullfile(srcdir,files{find(idx,1)});
    if sum(idx)>1; warning('%s daugther files found: used %s',num2str(sum(idx)),daughterfile); end

    % return to original dir
    cd(thisdir)
end