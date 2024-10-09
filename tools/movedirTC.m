function movedirTC(sourceFolder,destinationFolder)
    if exist(destinationFolder,'dir')
        warning(['Overwriting ',destinationFolder])
        rmdir(destinationFolder,"s")
    end
    
    if ispc
        system(sprintf('move "%s" "%s"', sourceFolder, destinationFolder));
    elseif isunix
        system(sprintf('mv "%s" "%s"', sourceFolder, destinationFolder));
    else
        error('Unsupported operating system');
    end
end