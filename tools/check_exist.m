function check_exist(FileIn)
    if ~exist("FileIn")
        str = ['Not found: ', FileIn];
        error(str)
    end
end