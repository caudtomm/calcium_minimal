function movie = add_badperiods_as_nans(movie,badperiods,meta)
% add nans in the place of bad periods

if isempty(badperiods);return;end

tmp = nan(meta.height,meta.width,meta.numberframes);
bptab = [0 0 0; badperiods; 0 meta.numberframes+1 0];
for i_bp = 1:size(badperiods,1)+1
    idx = sum(bptab(1:i_bp+1,2)) - sum(bptab(1:i_bp,3));
    dur = bptab(i_bp+1,2)-bptab(i_bp,3)-1;
    tmp(:,:,bptab(i_bp,3)+1:bptab(i_bp+1,2)-1) = ...
        movie(:,:,idx-dur-i_bp+1:idx-i_bp);
end
movie = tmp; clear tmp