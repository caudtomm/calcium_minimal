function movedirTC(sourceFolder,destinationFolder)
    
    if ispc
        if exist(destinationFolder,'dir')
            warning(['Overwriting ',destinationFolder])
            destinationFolder_tmp = strcat(destinationFolder,'_tmp');
            system(sprintf('move "%s" "%s"', sourceFolder, destinationFolder_tmp));
            rmdir(destinationFolder,"s")
            system(sprintf('move "%s" "%s"', destinationFolder_tmp, destinationFolder));
        else
            system(sprintf('move "%s" "%s"', sourceFolder, destinationFolder));
        end
    elseif isunix
        if exist(destinationFolder,'dir')
            warning(['Overwriting ',destinationFolder])
            destinationFolder_tmp = strcat(destinationFolder,'_tmp');
            system(sprintf('mv "%s" "%s"', sourceFolder, destinationFolder_tmp));
            rmdir(destinationFolder,"s")
            system(sprintf('mv "%s" "%s"', destinationFolder_tmp, destinationFolder));
        else
            system(sprintf('mv "%s" "%s"', sourceFolder, destinationFolder));
        end
    else
        error('Unsupported operating system');
    end
end