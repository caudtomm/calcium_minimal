function [putative_badperiods, putative_isbadframe] = ...
    predict_badperiods(data_raw, template, xcorr_th, activity_th, pl)
arguments
    data_raw Movie
    template double = data_raw.timeavg
    xcorr_th double = .5
    activity_th double = .95 % quantile
    pl logical = false
end

% initialize
movie = data_raw.stack;
framecount = data_raw.nfr;
xcorrvals = nan(framecount,1);
putative_isbadframe = false(framecount, 1);
putative_badperiods = [];

% take out high values (probably activity)
movie(movie > quantile(movie(:), activity_th)) = nan;
tic
for i_f = 1:framecount
    if mod(i_f,50)==0; disp(num2str(i_f)); end
    movie(:,:,i_f) = fillmissing2(movie(:,:,i_f),"linear");
end
toc

% first prediction based on norm 2d cross-corr
tic
for i_f = 1:framecount
    if mod(i_f,50)==0; disp(num2str(i_f)); end
    xcorrvals(i_f) = mean(normxcorr2(template, movie(:,:,i_f)),"all");
    putative_isbadframe(i_f) = xcorrvals(i_f) <= xcorr_th;
end
toc

if pl
    y = [squeeze(xcorrvals), squeeze(double(putative_isbadframe))];

    figure
    plot(y)
end


end
