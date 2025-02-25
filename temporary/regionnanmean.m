function S = regionnanmean(data, L)
% S = regionnanmean(data, L)
% data: M-by-N matrix
% L: integer indicating region size. size(data,1) must be a multiple of L.
%
% Takes the nanmean of L-long consecutive regions along dim-1.
%
% S is a M/L-by-N matrix
%

nregions = size(data,1)/L;
mask = repelem([1:nregions]',L);

data3d = [];
for i_reg = 1:nregions; data3d(:,:,i_reg) = data(mask==i_reg, :); end

S = shiftdim(nanmean(data3d,1))';




end