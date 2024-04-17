function [movie_without_bp, adapted_badperiods] = subtract_badperiods(movie,badperiods)
% Subtract bad periods

numberframes = size(movie,3);

% adapt bad periods to a sub-trial movie
adapted_badperiods = [];
for i_bp = 1:size(badperiods,1)
    thisbp = badperiods(i_bp,:);
    if thisbp(2)>numberframes; continue; end
    if thisbp(3)>numberframes; thisbp(3)=numberframes; end
    adapted_badperiods = [adapted_badperiods; thisbp];
end

for i_bp = 1:size(badperiods,1)
    thisbp = adapted_badperiods(size(adapted_badperiods,1)-i_bp+1,:);
    movie(:,:,thisbp(2):thisbp(3)) = [];
end
movie_without_bp = movie;