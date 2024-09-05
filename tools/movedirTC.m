function movedirTC(sourceFolder,destinationFolder)
    if exist(destinationFolder,'dir')
        warning(['Overwriting ',destinationFolder])
        rmdir(destinationFolder,"s")
    end
    eval(sprintf('!move "%s" "%s"',sourceFolder,destinationFolder))
end