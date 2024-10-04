function cdtol(path)
% tolerant cd, works with Unix and Windows

cd(fullfile(strrep(path,'\','/')))

end
