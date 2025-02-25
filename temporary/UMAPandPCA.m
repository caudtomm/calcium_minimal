function [a, reduction2, X, stims] = UMAPandPCA(experiment, todo_fish, stims2use, trialn2usein,method)

%% initialize global vars
global ntrials L ncells labs s C pathout tframe fs trialn2use stims2use X

% wish to use pca scores as input to UMAP?
usepca = false;
% time period to include (relative to [stim_on, stim_off]) [s]
tframe = [1 0];

[a, X, stims, traces,reltrialn] = selectData(experiment, todo_fish,stims2use,trialn2usein,method);
switch method
    case 'time'
        a = nanzscore(a,[],2);
    case 'trialavg'
        a = nanzscore(a,[],2);
end

%% pca
[coef,score,latent,tsquared,explained,mu] = pca(a);
figure; subplot(121); imagesc(score), title('score'); subplot(122); imagesc(coef), title('coef');

%% characterize PC space

% how much of the variance does each PC explain?
nshuffle = 50;
shuffleexplained = [];
for i_sh = 1:nshuffle
    shufflecoef = coef(randperm(size(coef,1)),:);
    shufflescore = a*shufflecoef;
    shufflepcvar = std(shufflescore,[],1,'omitnan').^2;
    shuffleexplained(i_sh,:) = (shufflepcvar*100)./sum(shufflepcvar);
end
shufflemean = nanmean(shuffleexplained,1);
shufflestd = std(shuffleexplained);
pcthr = movmean(shufflemean + shufflestd*2,5)';
maxPC = find(explained<=pcthr,1)-1

figure; subplot(121); imagesc(shufflescore), title('shufflescore');
subplot(122); imagesc(shufflecoef), title('shufflecoef');

figure; hold on
y = shufflemean; t = 1:length(y);
curve1 = y + shufflestd; curve2 = y - shufflestd;
h = patch([t,fliplr(t)],[curve1, fliplr(curve2)],'g','FaceAlpha',.3,'EdgeColor','none');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot(t,y,'Color','g','LineWidth',2)
plot(t,explained,'Color','b','LineWidth',2)
plot(t,pcthr,'Color','r','LineStyle','--','LineWidth',1)
xlim([t(1),t(end)])
set(gcf, 'color', 'none');    
set(gca, 'color', 'none', 'XColor','w', 'YColor','w', 'ZColor','w');
legend({'shuffle','data','thres (2 STD)'},'Box','on','color','none','Location','best','EdgeColor','w','TextColor','w')
ylabel('explained variance')
xlabel('PC #')


%% visualize PC space


switch method
    case 'time'
        score2 = permute(traceFormat(score,L),[1 3 2]);
        fig3dplot(score2(:,1:3,:),{'PCA123','PC#1','PC#2','PC#3'},[],0);
        % fig3dplot(score2(:,2:4,:),{'PCA234','PC#2','PC#3','PC#4'},[],0);
        
        figscatter(score2(:,1:2,:),{'PCA12','PC1','PC2',''},1)
        figscatter(score2(:,2:3,:),{'PCA23','PC2','PC3',''},1)
        figscatter(score2(:,[1,3],:),{'PCA13','PC1','PC3',''},1)
        
        
        fig3dplot(score2(:,[1,3,4],:),{'PCA134','PC#1','PC#3','PC#4'},[],0);
    case 'trialavg'

        % or
        
        h = figure
        for i_trial = 1:size(score,1)
        hp = scatter(score(i_trial,1),score(i_trial,2),[],C(X(i_trial),:),'filled');
        hold on
        end
        
        % cosmetics / labels
        set(gcf, 'color', 'none');    
        set(gca, 'color', 'none', 'XColor','w', 'YColor','w', 'ZColor','w');
        xlabel('PC1')
        ylabel('PC2')
        grid off
        axis square
        set(h,'Position',[573 438 300 250])
        
        
        % kmeans clustering
        metric = 'cosine'
        clust = kmeans(score-mean(score),4,'distance',metric)+1; % +1 to avoid black on black markers
        h = figure
        for i_trial = 1:size(score,1)
        hp = scatter(score(i_trial,1),score(i_trial,2),[],C(clust(i_trial),:),'filled');
        hold on
        end
        
        % cosmetics / labels
        set(gcf, 'color', 'none');    
        set(gca, 'color', 'none', 'XColor','w', 'YColor','w', 'ZColor','w');
        xlabel('PC1')
        ylabel('PC2')
        grid off
        axis square
        set(h,'Position',[573 438 300 250])
        
        % distance matrix used by kmeans
        x = squareform(pdist(score,metric));
        figure; imagesc(1-x)
        yticks(1:numel(X))
        xticks(1:numel(X))
        yticklabels(stims(X))
        xticklabels(stims(X))
        axis square
        colormap('hot')
        set(gcf, 'color', 'none');    
        set(gca, 'color', 'none', 'XColor','w', 'YColor','w', 'ZColor','w');
        colorbar('Color','w')

end

% store clustering accuracy to array (numel: fish)
% ACC = [ACC;cluster_acc(X,clust)]

%% print clustering accuracy over fish

% y = ACC;

scalefact = .3;

bias = cumsum(ones(size(y)),2)-1;
x = bias + scalefact*(rand(size(y))-.5);
h = figure;
chancelev = 1/numel(unique(X)); % assuming each repeats an equal num of times
hold on
line([min(x,[],'all')-.1,max(x,[],'all')+.1],repelem(chancelev,2), ...
    'Color','w','LineStyle',':','LineWidth',2)
for i = 1:size(y,2)
    line([min(x(:,i))-.2,max(x(:,i))]+.2,repelem(mean(y(:,i)),2), ...
        'Color','r','LineStyle',':','LineWidth',2)
end
for i = 1:size(y,1)
line(x(i,:),y(i,:),'Color','w','LineStyle','-','LineWidth',.5)
end
scatter(x,y,'r','filled')
set(gcf, 'color', 'none');
set(gca, 'color', 'none', 'XColor','w', 'YColor','w', 'ZColor','w', 'FontSize', 12);
xticks([])
yticks(linspace(0,1,6))
xlim([min(x,[],'all')-.3,max(x,[],'all')+.3])
ylim([0,1])
ylabel('kmeans accuracy')
set(h,'Position',[573 438 50+50*size(y,2) 450])


%% umap 
if usepca
    new_score = score(:,1:min([size(score,2), maxPC]));
else
    new_score = a;
end

% temporarily remove nan rows to be able to run umap (will be reinserted)
idx = isnan(mean(new_score,2));
new_score(idx,:) = [];

[reduction, umap, clusterIdentifiers, extras]=run_umap(new_score, ...
    'metric','euclidean', ...
    'min_dist',.25, ... % .8 for cosine, .25 for euclidean
    'n_components',2, ...
    'n_neighbors', 199, ...
    'init','spectral');
close

% reinsert nan rows
new_score = insert_rows(new_score,find(idx));
reduction = insert_rows(reduction,find(idx));

switch method
    case 'time'
        reduction2 = nan(size(traces,1),size(reduction,2),length(X));
        for i_trial = 1:size(traces,3)
            idx = (i_trial-1) * size(traces,1) + [1:size(traces,1)];
            reduction2(:,:,i_trial) = reduction(idx,:);
        end
    case 'trialavg'
        reduction2 = [];
        for i=1:size(reduction,1)
            reduction2(:,:,i)=reduction(i,:);
        end
end

% visualize UMAP space
figscatter(reduction2,{'UMAP','','',''},0)




%%
% line plot (scale by trel trial num)
fig3dplot(reduction2,{'UMAP_thickness','UMAP#1','UMAP#2','UMAP#3'}, ...
    numel(trialn2use)./(reltrialn*2),0);

% line plot selected cell
cellid = randi(size(traces,2),1);
thisisactive = squeeze(traces(:,cellid,:))>=.2;
[~,FileOut] = fig3dplot(reduction2,{'UMAP','UMAP#1','UMAP#2','UMAP#3'},[],0);
for i_trial = 1:ntrials
    hsc = scatter3(reduction(idx,1),reduction(idx,2),reduction(idx,3), ...
        40,'w','filled','DisplayName',['cell #' num2str(cellid)]);
    if i_trial > 1; hsc.HandleVisibility = 'off'; end
end
animate3D(FileOut,3);

% line plot smooth
fig3dplot(movmean(reduction2,3,1),{'UMAP_smooth','UMAP#1','UMAP#2','UMAP#3'},[],0);

% 
% % dUMAP
% reddiff = diff(movmean(reduction2,5,1),[],1);
% fig3curves(reddiff,{'dUMAP_smooth','dUMAP#1','dUMAP#2','dUMAP#3'});
% 
% % zeroed-out UMAP
% reddiffitg = cumsum(reddiff,1);
% fig3curves(reddiffitg,{'relUMAP_smooth_curves','relUMAP#1','relUMAP#2','relUMAP#3'});
% fig3curves(reddiffitg,{'relUMAP_smooth_curvesavg','relUMAP#1','relUMAP#2','relUMAP#3'},1);
% fig3dplot(reddiffitg,{'relUMAP_smooth_3d','relUMAP#1','relUMAP#2','relUMAP#3'},[],1);
% 
% y = permute(reddiffitg,[1,3,2]);
% d = nan(size(y,2),size(y,2),size(y,1));
% for i_trial1 = 1:size(y,2)
%     reftrial = squeeze(y(:,i_trial1,:));
%     for i_trial2 = 1:size(y,2)
%         thistrial = squeeze(y(:,i_trial2,:));
%         d(i_trial1,i_trial2,:) = distD(reftrial',thistrial');
%     end
% end
% c = nanmean(d,3);
% [~,idx] = sort(X);
% h = figure; imagesc(1-c(idx,idx)); axis square
% FileOut = fullfile(pathout,'relUMAPDistAVG.fig');
% savefig(h,FileOut)
% 
% 
% % rotate and scale to unit vector
% y = movmean(reduction2(1:end,:,:),5,1);
% rs = nan(size(y));
% trialref = [];
% for i_trial = 1:size(y,3)
%     thistrial = squeeze(y(:,:,i_trial));
%     rs(:,:,i_trial) = rotatescale3(thistrial,trialref);
%     if i_trial == 1; trialref = squeeze(rs(:,:,i_trial)); end
% end
% rs_smooth = rs;
% figlabs = {'rotnormUMAP_smooth','UMAP#1','UMAP#2','UMAP#3'};
% [h,FileOut] = fig3dplot(rs_smooth,figlabs, ...
%     numel(trialn2use)./(reltrialn*2),0);
% scatter3([0,1],[0,0],[0,0],50,'w','filled','HandleVisibility','off');
% FileOut2 = fullfile(pathout,strcat(figlabs{1},'.fig'));
% if s; savefig(h,FileOut2); end
% animate3D(FileOut,3);
% fig3curves(rs_smooth,{'rotnormUMAP_smooth_curves','rotnormUMAP#1','rotnormUMAP#2','rotnormUMAP#3'},1);
% 
% % quantify relative distances
% y = permute(rs,[1,3,2]);
% d = nan(size(y,2),size(y,2),size(y,1));
% for i_trial1 = 1:size(y,2)
%     reftrial = squeeze(y(:,i_trial1,:));
%     for i_trial2 = 1:size(y,2)
%         thistrial = squeeze(y(:,i_trial2,:));
%         d(i_trial1,i_trial2,:) = distD(reftrial',thistrial');
%         %d(i_trial1,i_trial2,:) = mean(convn(thistrial, reftrial(end:-1:1,end:-1:1),2));
%     end
% end
% 
% c = nanmean(d,3);
% [~,idx] = sort(X);
% h = figure; imagesc(1-c(idx,idx)); axis square
% FileOut = fullfile(pathout,'rotnormDistAVG.fig');
% savefig(h,FileOut)




