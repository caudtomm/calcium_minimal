function structout = getFileNameSpecs(filename)
filename = char(filename);

structout.filename = filename;
[structout.fname, structout.orig_fext, structout.orig_fpath] = getFName(filename);

% just out of laziness
fname = structout.fname;

structout.owner = extractval(fname); fname = extractAfter(fname,'_');
structout.date = extractval(fname); fname = extractAfter(fname,'_');
structout.subject_line = extractval(fname); fname = extractAfter(fname,'_');
structout.subject = extractval(fname); fname = extractAfter(fname,'_');
structout.region = extractval(fname); fname = extractAfter(fname,'_');
structout.stim_type = extractval(fname); fname = extractAfter(fname,'_');
structout.method = extractval(fname); fname = extractAfter(fname,'_');

if ~isempty(fname) % correct for double '_'
    while true
        if strcmp(fname(1),'_')
            fname = fname(2:end);
        else
            break
        end
    end
end

structout.trial_num = extractval(fname,'num');
if ~isempty(structout.trial_num); fname = extractAfter(fname,'_'); end

structout.extra = fname;

end

%%

function val = extractval(inchar, outtype)

if ~exist("outtype","var"); outtype = ''; end

switch outtype
    case 'num'
        if isempty(extractBefore(inchar,'_'))
            val = str2num(inchar);
        else
            val = str2num(extractBefore(inchar,'_'));
        end
    otherwise
        if isempty(extractBefore(inchar,'_'))
            val = inchar;
        else
            val = extractBefore(inchar,'_');
        end
end

end