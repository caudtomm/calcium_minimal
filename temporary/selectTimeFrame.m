function M = selectTimeFrame(traces, interval, L)
% this funtion is really thrown together in a rush, sorry...
% 
% returns a matrix M with n dim = n dim of input traces.
% dim 1 is sized down to include only values found within integer multiples
% of the input timeframe length 'L', plus the input 'interval', whenever
% entirely contained within dim 1.
% works for 1,2,3 dim traces.
% 
% input
% traces : dim 1 is time
% interval : [start_idx : end_idx]
% (L) : timeframe length (ex. length of one trial in frames) - default : size(traces,1) 
% 
% M = selectTimeFrame(traces, period, L)
% 
% ex.
% M = selectTimeFrame(magic(10), 2:3, 5);
% 
% Version:
% 16 Apr 2020 - Tommaso Caudullo
% 

switch nargin 
    case {0, 1}
        error('Not enough input arguments!')
    case 2
        L = size(traces,1);
    case 3
    otherwise
        error('Too many input arguments!')
end

% initialize
M = traces;

% execute
reps = floor(size(traces,1)/L);
mask = L * repelem(0:reps-1,length(interval)) + repmat(interval,[1,reps]);

M = traces(mask,:,:);





end