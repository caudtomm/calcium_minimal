%%
baseline_tag = {'noodor', 'baseline', 'spont.'};
this_group = {'previousnaive';'naive';'trained1';'trained2';'trained1alt';'uncoupled'};
cat.traces = [];
todo_fish = [];
for i_fish = 1:numel(experiment.series)
    if ismember(experiment.series{i_fish}.group,this_group)
        todo_fish = [todo_fish;i_fish];
    end
end

for i_fish = 1:numel(todo_fish)
    %%
    data = experiment.series{i_fish}.data;
    tmp = data.tracesdesdn;
    fs = data.meta.framerate;
%     ds = data.meta.downsample;
    ds = 1;


    
%     FindTopUnits
%     tmp = selectCells(tmp,data.L/ds,alltopunits');
    tmp = traceFormat(tmp,data.L/ds);
    

    % zscore each cell, based of its mean and variance throughout the
    % experiment (all trials)
    tmp2 = [];
    for i=1:size(tmp,3)
    tmp2 = [tmp2;tmp(:,:,i)];
    end
    for i = 1:size(tmp2,2)
    tmp2(isnan(tmp2(:,i)),i) = nanmean(tmp2(:,i));
    end
    tmp2 = zscore(tmp2);
    tmp3 = [];
    for i=1:size(tmp,3)
    tmp3(:,:,i) = tmp2(1300*(i-1)+[1:1300],:);
    end
    tmp = tmp3;

    % component to isolate
    tmp(tmp<0) = nan;



    % take out baseline trials
    tmp(:,:,ismember(data.stim_type, baseline_tag)) = [];

    %zscore
%     T = size(tmp,1);
%     tmp = zscore(traceFormat(permute(tmp,[1,3,2])),[],[],'omitnan');
%     tmp = permute(traceFormat(tmp,T),[1,3,2]);
    
    % select odor window [-5, +20] seconds
    interval = [floor(nanmax([1, (data.stim_on_sec - 4)*fs/ds])) : ...
        floor(nanmin([size(tmp,1), (data.stim_off_sec + 20)*fs/ds]))];
%     interval = [floor(nanmax([1, (data.stim_off_sec + 94)*fs/ds])) : ...
%         floor(nanmin([size(tmp,1), (data.stim_off_sec + 118)*fs/ds]))];
%     interval = [floor(nanmax([1, (data.stim_on_sec - 4)*fs/ds])) : ...
%         floor(nanmin([size(tmp,1), (data.stim_off_sec + 104)*fs/ds]))];
    tmp = tmp(interval,:,:);

   
%     tmp = tmp(:,:,[1:5]+5*2);
    % take avg over cells
    tmp = squeeze(nanmean(tmp,2));
    
    % normalize
%     tmp = tmp - quantile(tmp,.03); %- nanmean(tmp(1:10,:,:),[1,2]);
    
    % take out bad trials
    
%     % plot
%     figure; plot(tmp)
    
    %% store to struct
    if size(cat.traces,1)>size(tmp,1)
        tmp = [tmp; nan(size(cat.traces,1)-size(tmp,1), size(tmp,2))];
    elseif size(cat.traces,1)<size(tmp,1)
        cat.traces = [cat.traces; nan(size(tmp,1)-size(cat.traces,1), size(cat.traces,2))];
    end

    for i = 1:size(tmp,2)
        tmp(isnan(tmp(:,i)),i) = nanmean(tmp(:,i));
%         tmp(:,i) = zscore(tmp(:,i));
    end
    
    cat.traces = [cat.traces,nanmean(tmp,2)]; %
end


% x = [32,34,37,64,132];
% for i = 1:length(x)
% cat.traces(:,x(end-i+1)) = [];
% end

avg = nanmean(cat.traces(1:floor(4*fs/ds),:),[1,2])
% avg = 0;
%% plot avg dFoverF
t = linspace(-4,40,size(cat.traces,1))
y = nanmean(cat.traces,2)-avg;
figure; hold on
curve1 = y + std(cat.traces,[],2,'omitnan'); curve2 = y - std(cat.traces,[],2,'omitnan');
h = patch([t,fliplr(t)],[curve1; fliplr(curve2')'],'b','FaceAlpha',.3,'EdgeColor','none')
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot(t,y, 'LineWidth', 5, 'Color', 'b')
% line([0,0],[-.01,.015],'Color','r','LineWidth',2,'LineStyle','--')
% line([20,20],[-.01,.015],'Color','r','LineWidth',2,'LineStyle','--')

xlabel('time from stimulus onset [s]')
ylabel('norm. intensity change')
axis tight

set(gcf, 'Position', [50 50 400 280]);
set(gca, 'color', 'none', 'XColor','w', 'YColor','w', 'ZColor','w');
set(gcf, 'color', 'none'); 

l = legend('baseline','odor','pos.','neg.'); legend('boxoff')
l.TextColor = 'w';


%% across trials

ytmp = squeeze(nanmean(cat.traces,1));

numtrials_fish = [];
for i_series = 1:numel(todo_fish)
    thisfish = todo_fish(i_series);
    data = experiment.series{thisfish}.data;
    numtrials_fish = [numtrials_fish;
                      numel(data.trial_num(~ismember(data.stim_type, baseline_tag)))];
end 
[m,i] = max(numtrials_fish);
y = nan(m,numel(todo_fish));
data = experiment.series{i}.data;
for i_series = 1:numel(todo_fish)
    fish_range = sum(numtrials_fish(1:i_series-1))+[1:numtrials_fish(i_series)];
    if isempty(fish_range); fish_range = 1:numtrials_fish(i_series); end
    if isempty(fish_range); continue; end
    y(1:numtrials_fish(i_series),i_series) = ytmp(fish_range);
end
ymean = nanmean(y,2);
y = y./nanmean(ymean);
ymean = ymean./nanmean(ymean);
t = [1:(length(ymean))]';

h = figure; set(gcf, 'Position', [50 50 1200 500],'color','k')
subplot(121); scatter(t,y,'w','filled'); hold on
% for i=1:size(y,2); plot(t,y(:,i),'w'); end
scatter(t,ymean,100,'r','filled')
xticks(1:max(t)); xticklabels(data.stim_type(2:34))
ylabel('norm. intensity','Color','w','FontSize',16);xlim([0,max(t)+1])
set(gca, 'color', 'k', 'XColor','w', 'YColor','w', 'ZColor','w');
idx = data.idx_by_stim_type;
idx(ismember(data.stim_type(idx),baseline_tag))=[];
idx(ismember(data.stim_type(idx),{'ACSF'}))=[];
subplot(122); scatter(t,y(idx-1,:),'w','filled'); hold on
% for i=1:size(y,2); plot(t,y(idx,i),'w'); end
scatter(t,ymean(idx-1),100,'r','filled')
xticks(1:max(t)); xticklabels(data.stim_type(idx));
ylabel('norm. intensity','Color','w','FontSize',16); xlim([0,max(t)+1])
set(gca, 'color', 'k', 'XColor','w', 'YColor','w', 'ZColor','w');



h = figure; set(gcf, 'Position', [50 50 600 500],'color','k')
scatter(t,y,'w','filled'); hold on
% for i=1:size(y,2); plot(t,y(:,i),'w'); end
scatter(t,ymean,100,'r','filled')
xticks(1:max(t)); xticklabels(data.stim_type(2:34))
ylabel('norm. intensity','Color','w','FontSize',16);xlim([0,max(t)+1])
set(gca, 'color', 'k', 'XColor','w', 'YColor','w', 'ZColor','w');

kruskalwallis(y(1:5,:)')

h = figure; set(gcf, 'Position', [50 50 600 500],'color','k')
idx = data.idx_by_stim_type(1:32)-1;
scatter(t,y(idx,:),'w','filled'); hold on
% for i=1:size(y,2); plot(t,y(idx,i),'w'); end
scatter(t,ymean(idx),100,'r','filled')
xticks(1:max(t)); xticklabels(data.stim_type(idx+1));
ylabel('norm. intensity','Color','w','FontSize',16); xlim([0,max(t)+1])
set(gca, 'color', 'k', 'XColor','w', 'YColor','w', 'ZColor','w');



%% across blocks
baseline_tag = {'noodor', 'baseline', 'spont.'};

% blocktrials_idx = {[1:5]+5*0; ... % block 1
%                    [1:5]+5*1; ... % block 2
%                    [1:5]+5*2; ... % block 3
%                    [1:7]+5*3; ... % block TRP
%                    [1:5]+5*4+2; ... % block 4
%                    [1:5]+5*5+2}; % block 5
blocktrials_idx = {[1:5]+6*0; ... % block 1
                   [1:5]+6*1; ... % block 2
                   [1:5]+6*2; ... % block 3
                   [1:5]+6*3; ... % block 4
                   [1:5]+6*4; ... % block 5
                   [1:5]+6*5}; % block 6

i_block = 6;


cat.traces = [];

for i_fish = 1:numel(experiment.series)
    %%
    data = experiment.series{i_fish}.data;
    tmp = data.tracesdn;
    fs = data.meta.framerate;
%     ds = data.meta.downsample;
    ds = 1;
    tmp(tmp<.2) = .2;
    
%     FindTopUnits
%     tmp = selectCells(tmp,data.L/ds,alltopunits');
    tmp = traceFormat(tmp,data.L/ds);
    
    % take out baseline trials
    tmp(:,:,ismember(data.stim_type, baseline_tag)) = [];
    
    % select odor window [-5, +20] seconds
    interval = [floor(nanmax([1, (data.stim_on_sec - 4)*fs/ds])) : ...
        floor(nanmin([size(tmp,1), (data.stim_off_sec + 20)*fs/ds]))];
%     interval = [floor(nanmax([1, (data.stim_off_sec + 30)*fs/ds])) : ...
%         floor(nanmin([size(tmp,1), (data.stim_off_sec + 54)*fs/ds]))];
    tmp = tmp(interval,:,:);
   
    tmp = tmp(:,:,blocktrials_idx{i_block});
    % take avg over cells
    tmp = squeeze(nanmean(tmp,2));
    
    % normalize
    tmp = tmp - quantile(tmp,.03); %- nanmean(tmp(1:10,:,:),[1,2]);
    
    % take out bad trials

%     % plot
%     figure; plot(tmp)
    
    %% store to struct
    if size(cat.traces,1)>size(tmp,1)
        tmp = [tmp; nan(size(cat.traces,1)-size(tmp,1), size(tmp,2))];
    elseif size(cat.traces,1)<size(tmp,1)
        cat.traces = [cat.traces; nan(size(tmp,1)-size(cat.traces,1), size(cat.traces,2))];
    end

    for i = 1:size(tmp,2)
        tmp(isnan(tmp(:,i)),i) = nanmean(tmp(:,i));
        tmp(:,i) = zscore(tmp(:,i));
    end

    cat.traces = [cat.traces, tmp];
end

% x = [5,7,10,24];
% for i = 1:length(x)
% cat.traces(:,x(end-i+1)) = [];
% end

avg = nanmean(cat.traces(1:floor(4*fs/ds),:),[1,2])

%% plot avg dFoverF
figure; hold on
t = linspace(-4,40,size(cat.traces,1))
bias = (numel(blocktrials_idx)-i_block)*.03;
y = nanmean(cat.traces,2)-avg + bias;
curve1 = y + std(cat.traces,[],2,'omitnan'); curve2 = y - std(cat.traces,[],2,'omitnan');
h = patch([t,fliplr(t)],[curve1; fliplr(curve2')'],'b','FaceAlpha',.3,'EdgeColor','none')
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot(t,y, 'LineWidth', 5, 'Color', 'b')
line([0,0],[-.01,.2],'Color','r','LineWidth',2,'LineStyle','--')
line([20,20],[-.01,.2],'Color','r','LineWidth',2,'LineStyle','--')

xlabel('time from stimulus onset [s]')
%ylabel('mean dF/F')
axis tight

set(gcf, 'Position', [50 50 400 800]);
set(gca,'ytick',[])
set(gca, 'color', 'none', 'XColor','w', 'YColor','w', 'ZColor','w');
set(gcf, 'color', 'none'); 

%% 













