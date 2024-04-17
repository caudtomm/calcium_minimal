function [movie_subtracted, baseline] = subtract_movie_baseline(movie,framerate)
% baseline defined as the 3% quantile of all pixel values per frame (50 sec rolling avg)
baseline = movmean(quantile(movie,0.03,[1,2]),floor(50*framerate),'omitnan'); 

baselineMAT = repmat(baseline,[size(movie,1),size(movie,2)]);
movie_subtracted = movie - baselineMAT;

baseline = squeeze(baseline);
