
%% plot trial distance averaged across fish

% todo_fish = [6:7,11:12,14,17,18,21:23];
todo_fish = find(~strcmp(experiment.summaryTable.group,'naïve'));
% todo_fish = 20;
% todo_fish = 1:17

sec_range = [1,20]; % [+stim_on, +stim on]
crange = [0.6 1];
% blockstructure = [5;5;5;7;5;5];
blockstructure = [6;6;6;6;6;3];

data = experiment.series{todo_fish(1)}.data;
% sortorder = 1:numel(data.stim_type);
sortorder = data.idx_by_stim_type;
labs = data.stim_type(sortorder);

C = nan(numel(sortorder),numel(sortorder),numel(todo_fish));
for i_fish = 1:numel(todo_fish)
    data = experiment.series{todo_fish(i_fish)}.data;
    fs = data.meta.framerate;
    this_sortorder = [1:numel(data.trials)]';
    interval = [floor((data.stim_on_sec+sec_range(1))*fs) : ...
              1+floor((data.stim_on_sec+sec_range(2))*fs)];
    c = distTraceData(data,this_sortorder,data.L,[],'cosine',[1:data.N]',interval,0,0);
    C(1:size(c,1),1:size(c,2),i_fish) = c;
%     figure; imagesc(C(1:16,1:16,i_fish)), title(num2str(i_fish))
end

% plot full avg intertrial distance matrix
c = nanmean(C(sortorder,sortorder,:),3);
figure; imagesc(c)
axis square; hold on
xticks(1:size(c,1)); xticklabels(labs); xtickangle(90)
yticks(1:size(c,1)); yticklabels(labs)
xlabel('Stimulus type')
ylabel('Stimulus type')
title(['Similarity: cosine ', num2str(interval(1)/fs-data.stim_on_sec), '-', ...
    num2str(interval(end)/fs-data.stim_on_sec), ' s'], 'Color','w')
caxis([quantile(c(:),.001), quantile(c(:),.999)]);
caxis(crange)
colormap(['hot'])
colorbar('Color','w')
set(gca, 'color', 'none', 'XColor','w', 'YColor','w', 'ZColor','w');
set(gcf, 'color', 'none'); 
% set(gcf, 'Position', [50 50 200 400]);
hold off

%% plot average distance between trials of the same odor for each block

ref = MakeRef34MAT(experiment.series{todo_fish(1)}.data.stim_type,0);

% select only trials of the same odors, but exclude diagonal elements and
% spontaneous activity trials
mask = ref; mask(mask==0)=nan;
c = C.*mask; %c = c + diag(nan(34,1) - diag(c));
c = c(2:end-1,2:end-1);

% take mean of each block
out = nan(numel(blockstructure));
for row = 1:numel(blockstructure)
    for col = 1:numel(blockstructure)
        rowrange = sum(blockstructure(1:row-1)) + [1:blockstructure(row)];
        colrange = sum(blockstructure(1:col-1)) + [1:blockstructure(col)];
        out(row,col) = nanmean(c(rowrange,colrange),'all');
    end
end

same_odors = out; same_odors(isnan(same_odors)) = 0;
c = out;
figure; imagesc(c)
axis square; hold on
xticks(1:size(c,1)); 
yticks(1:size(c,1)); 
xlabel('Block #')
ylabel('Block #')
title('Same odor similarity', 'Color','w')
caxis(crange);
colormap(['hot'])
colorbar('Color','w')
set(gca, 'color', 'none', 'XColor','w', 'YColor','w', 'ZColor','w','FontSize',14);
set(gcf, 'color', 'none'); 
set(gcf, 'Position', [50 50 350 400]);
hold off



%% plot average distance between trials of the different odors for each block

ref = MakeRef34MAT(experiment.series{todo_fish(1)}.data.stim_type,0);

% select only trials of the same odors, but exclude diagonal elements and
% spontaneous activity trials
mask = double(~logical(ref)); mask(mask==0)=nan;
c = C.*mask; %c = c + diag(nan(34,1) - diag(c));
c = c(2:end-1,2:end-1);

% take mean of each block
out = nan(numel(blockstructure));
for row = 1:numel(blockstructure)
    for col = 1:numel(blockstructure)
        rowrange = sum(blockstructure(1:row-1)) + [1:blockstructure(row)];
        colrange = sum(blockstructure(1:col-1)) + [1:blockstructure(col)];
        out(row,col) = nanmean(c(rowrange,colrange),'all');
    end
end

different_odors = out; different_odors(isnan(different_odors)) = 0;
c = out;
figure; imagesc(c)
axis square; hold on
xticks(1:size(c,1)); 
yticks(1:size(c,1)); 
xlabel('Block #')
ylabel('Block #')
title('Different odors similarity', 'Color','w')
caxis(crange);
colormap(['hot'])
colorbar('Color','w')
set(gca, 'color', 'none', 'XColor','w', 'YColor','w', 'ZColor','w','FontSize',14);
set(gcf, 'color', 'none'); 
set(gcf, 'Position', [50 50 350 400]);
hold off




%% plot average distance between trials of the same - different odors for each block

c = same_odors - different_odors;
figure; imagesc(c)
axis square; hold on
xticks(1:size(c,1)); 
yticks(1:size(c,1)); 
xlabel('Block #')
ylabel('Block #')
title('Same-different odors similarity', 'Color','w')
caxis([0 1]);
colormap(['hot'])
colorbar('Color','w')
set(gca, 'color', 'none', 'XColor','w', 'YColor','w', 'ZColor','w','FontSize',14);
set(gcf, 'color', 'none'); 
set(gcf, 'Position', [50 50 350 400]);
hold off

%% plot similarity of adjacent blocks (same-different diag+1 elements)

c = same_odors - different_odors;

figure; plot(diag(c,1),'LineWidth',2,'Color','r')

ylim([0 1])
xticks(1:size(c,1)); xticklabels({'1-2','2-3','3-4','4-5','5-6'}); xtickangle(45)
xlabel('Block #')
ylabel('same-different similarity')
title('Adjacent blocks', 'Color','w')
set(gca, 'color', 'none', 'XColor','w', 'YColor','w', 'ZColor','w','FontSize',14);
set(gcf, 'color', 'none'); 
set(gcf, 'Position', [50 50 350 400]);
box off

