function movie = add_badperiods_as_nans(movie_raw)
% add nans in the place of bad periods
arguments
    movie_raw Movie
end

% initialize
movie = movie_raw.stack;
badperiods = movie_raw.badperiods;
h = movie.h;
w = movie.w;
nfr = movie.nfr;

if isempty(badperiods);return;end

tmp = nan(h,w,nfr);
bptab = [0 0 0; badperiods; 0 nfr+1 0];
for i_bp = 1:size(badperiods,1)+1
    idx = sum(bptab(1:i_bp+1,2)) - sum(bptab(1:i_bp,3));
    dur = bptab(i_bp+1,2)-bptab(i_bp,3)-1;
    tmp(:,:,bptab(i_bp,3)+1:bptab(i_bp+1,2)-1) = ...
        movie(:,:,idx-dur-i_bp+1:idx-i_bp);
end
movie = tmp; clear tmp