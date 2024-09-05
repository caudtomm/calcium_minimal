
function structout = getFileNameSpecs(filename)
filename = char(filename);

structout.filename = filename;
[structout.fname, structout.orig_fext, structout.orig_fpath] = getFName(filename);

% just out of laziness
fname = structout.fname;

structout.owner = extractBefore(fname,'_'); fname = extractAfter(fname,'_');
structout.date = extractBefore(fname,'_'); fname = extractAfter(fname,'_');
structout.subject_line = extractBefore(fname,'_'); fname = extractAfter(fname,'_');
structout.subject = extractBefore(fname,'_'); fname = extractAfter(fname,'_');
structout.region = extractBefore(fname,'_'); fname = extractAfter(fname,'_');
structout.stim_type = extractBefore(fname,'_'); fname = extractAfter(fname,'_');
structout.method = extractBefore(fname,'_'); fname = extractAfter(fname,'_');

if ~isempty(fname)
    while true
        if strcmp(fname(1),'_')
            fname = fname(2:end);
        else
            break
        end
    end
end

if isempty(extractBefore(fname,'_'))
    structout.trial_num = str2num(fname);
else
    structout.trial_num = str2num(extractBefore(fname,'_'));
end

end