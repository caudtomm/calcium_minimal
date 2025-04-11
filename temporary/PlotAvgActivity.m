
%% load baseline

if ~exist('baseline', 'var')
    plane = robust_io('load',fullfile(FileIn_path,'defROIs', ...
        strcat(baselinefname,'_defROIs.mat')), 'plane').plane;
    [~, ~, ~, ~, baseline] = getTransients(plane{1}.timetraces);
end


%% average responses across neurons, per trial

figure(242); hold on; title('Average activity of common units')
L = data.Ldns;
traces = selectCells(data.tracesdns,L,cells');
factor = nanmean(baseline) * .9;
fs = data.meta.framerate/data.meta.downsample;
x = 0:1/fs:L/fs-1/fs;

col = optimalcolors(10);
count = 2;

for i_trial = 1:size(data.traces, 2)
    
%     c = rand(1,3);
    c = col(count,:);
    intercept = factor * (i_trial-1);
    
    % activity for this trial
    thistraces = reshape(traces(:,i_trial), [L numel(cells)])' ...
        - nanmean(baseline); % rows are cells  
    
    % get top 20 most variant cells during stimulus window + 10s
    interval = [floor(data.stim_on_sec*fs) : (data.stim_off_sec+10)*fs];
    cellstd = std(thistraces(:,interval),[],2);
    [~,idx] = sort(cellstd); idx = idx(end-19:end);
    
    % plot trace
    figure(242);
    b(count) = plot(x,intercept+nanmean(thistraces(idx,:),1), 'Color', c, 'LineWidth', 1.5);
    
end


xlabel('time [s]')

yticks(factor * [0:numel(data.trial_num)- 1] - nanmean(baseline)); yticklabels(data.stim_type)
axis tight
ylim([min(ylim)-nanmean(baseline) size(data.traces,2)*factor])
axis square
set(gcf,'Position',[5 5 900 700])

figure(242)
line([data.stim_on_sec data.stim_on_sec],ylim,'Color', 'r', 'LineStyle', '--', 'LineWidth',1)
line([data.stim_off_sec data.stim_off_sec],ylim,'Color', 'r', 'LineStyle', '--', 'LineWidth',1)


% legend(b,{'vpDp','dpDp'})


