function b = prompt_overwrite(filename)
    b = true;
    if ~exist(filename,'file'); return; end
    str = ['Pre-existing file found: \n',char(filename),'\nOverwrite? [y/n]'];
    b = askBoolean(str);
    if ~b; warning('Aborted overwrite.'); end
end
