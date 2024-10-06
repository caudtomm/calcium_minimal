function cdtol(thispath)
% tolerant cd, works with Unix and Windows

cd(fullfile(strrep(thispath,'\','/')))

end
