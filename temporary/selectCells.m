function tracesout = selectCells(tracesin, L, idx)
% tracesout = selectCells(tracesin, L, idx)
% 
% accepts data.traces format (cells are sequential in array, 
%                             dim 1 : time
%                             dim 2 : trials
% 
% L is length of each cell's trace (int)
% size(tracesin,1) must be an integer multiple of L
% 
% idx are the indices of all cells to KEEP (1:n integer array)
% 
% returns same trace format as input, but containing only data from cells
% indexed in (idx)
% 


x = [1:L]'+L*(squeeze(idx)-1);
tracesout = tracesin(x(:),:);
