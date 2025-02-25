function idx = getLinearCorrespondanceMap(eqitems)
% implements the expression n+sum_[i=1->N-1](i-1) iteratively on ascending
% unique pairs of entry elements
%
% idx = getLinearCorrespondanceMap(eqitems)
%
% Useful to get the indices of elements from M(triu(M)), where M:[NxN],
% excluding the diagonal. The elements are identified as all
% correspondances between 'eqitems' items on the original square matrix,
% counted only once.
%
% input:
%   eqitems : [nx1] , integer items whose correspondance indices are to be
%                     extracted. Each 1 <= eqitem(i) <= N
%

eqitems = sort(eqitems);
idx = [];
for i_n = 1:numel(eqitems)-1
    n = eqitems(i_n);
    
    for i_N = i_n+1:numel(eqitems)
        N = eqitems(i_N);
        
        tmp = 0;
        for k = 1:N-1
            tmp = tmp + k-1;
        end
        
        idx = [idx, n + tmp];
    end
end




end