function c = distTraceData(varargin)
% produces pairwise distance matrix according to method
% c = distTraceData(data, [sortorder], [L], [labels], [method], [cells], [interval], [pl])
%
%



%% initialize vars

if nargin==0; error('Not enough input arguments!'); end

data = varargin{1};

if isstruct(data)
    mode = 'structure'; 
    fs = data.meta.framerate; 
    traces = data.traces;
    defaultLen = data.L;
    defaultLabs = data.stim_type;
    defaultCells = data.common_units;
elseif isnumeric(data)
    mode = 'numeric';
    fs = 7.6667; % default framerate from galvo-galvo scanner
    traces = data;
    defaultLen = 1300;
    defaultLabs = [];
    defaultCells = [1:(size(data,1)/defaultLen)]';
else
    error('Input argument type not recognized!');
end

defaultDoZscore = false;
defaultPlot = true;

Defaults = {data, [1:size(traces,2)], defaultLen, defaultLabs, ...
    'cosine', defaultCells, 1:defaultLen, defaultDoZscore, defaultPlot};

if nargin > numel(Defaults); error('Too many input arguments.'); end

idx = ~cellfun('isempty',varargin);
Defaults(idx) = varargin(idx);



i = 2;

sortorder = Defaults{i}; i=i+1;
L = Defaults{i}; i=i+1;
labs = Defaults{i}; i=i+1;
method = Defaults{i}; i=i+1;
cells = Defaults{i}; i=i+1;
interval = Defaults{i}; i=i+1;
dozscore = Defaults{i}; i=i+1;
pl = Defaults{i};

clear i


interval(interval==0) = [];
if isempty(interval); error('Selected interval is too short!'); end



%% calc distances

% traces = selectCells(traces,L,transpose(find(ismember(cells,data.common_units))));
if dozscore
%     traces = traceFormat(zscore(traceFormat(traces,L)));
    traces = zscore(traces);
%     win_all = zscore(win_all);
end
traces = selectCells(traces,L,cells');
win_all = selectTimeFrame(traces(:,sortorder), interval, L);

avg_actvect = regionnanmean(win_all, length(interval));
% avg_actvect = avg_actvect - nanmean(avg_actvect,2);
% avg_actvect = zscore(avg_actvect,[],2);
win_d = squareform(pdist(avg_actvect', method));

c = win_d;
c = 1-c;



%% optional plotting (not saved)

if pl
    
    
%     % matrix diagonal is nans (for plotting only)
%     idx = logical(diag(ones(1,size(c,1))));
%     c(idx) = nan;
    
    
    
    figure; imagesc(c); axis square; hold on
    
%     if strcmp(mode, 'structure')
%         xticks(1:size(c,1)); xticklabels(sortorder); xtickangle(90)
        xticks(1:size(c,1)); xticklabels(labs(sortorder)); xtickangle(90)
        yticks(1:size(c,1)); yticklabels(labs(sortorder))
        xlabel('Trial number')
        ylabel('Stimulus type')
%     end
    
    title(['Distance: ', method, ' ', num2str(interval(1)/fs-data.stim_on_sec), '-', ...
    num2str(interval(end)/fs-data.stim_on_sec), ' s'])
    
    caxis([quantile(c(:),.1), quantile(c(:),.999)]);
%     caxis([.85 1])
    colorbar
    hold off
    
    
end



end