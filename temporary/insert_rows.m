function mat_out = insert_rows(mat_in,original_idx,original_rows)
% original_idx of type double (not boolean indices)

ndim = size(mat_in,2);

if ~exist('original_idx',"var"); mat_out = mat_in; return;  end
if ~exist('original_rows',"var"); original_rows = nan(length(original_idx),ndim); end
if size(original_rows,2)~=ndim; error('number of columns is inconsistent'); end

mat_out = nan( size(mat_in,1)+size(original_rows,1) , ndim );

for i = 1:size(mat_out,1)
    thisrow=[];
    position_in_idxmat = find(original_idx==i);
    if isempty(position_in_idxmat)
        thisrow = mat_in(1,:);
        mat_in(1,:) = [];
    else
        thisrow = original_rows(position_in_idxmat,:);
    end
    mat_out(i,:) = thisrow;
end