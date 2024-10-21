function win_d = getDistinRegn(traces, L, interval, method)


% function to distances off data (to package)


%interval = (80+10):(240+10); % set default
win_all = selectTimeFrame(traces, interval, L);
avg_actvect = regionnanmean(win_all, length(interval));

win_d = squareform(pdist(avg_actvect', method));



end




