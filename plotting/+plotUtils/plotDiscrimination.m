
%% plotting
figure;set(gcf,'Color','w')
for i_set = 1:nsets
    subplot(2,nsets,i_set); imagesc(correctLab(:,:,i_set))
    title(['training set ',num2str(i_set)])
    ylabel('fish')
end
for i_set = 1:nsets
    subplot(2,nsets,i_set+nsets); imagesc(correctLabSH(:,:,i_set))
end
xticks(1:numel(X)); xticklabels(stims(X));


C = optimalcolors(nsets+1); C = C(2:end,:);

h = figure; hold on
t = 1:ntrials;
for i_set = 1:nsets
    tmp = correctLab(:,:,i_set);
    ymean = 100* nanmean(tmp,1);
    ystd = 100* squeeze(std(tmp,[],1,'omitnan'));

    curve1 = ymean + ystd; curve2 = ymean - ystd;
    hp = patch([t,fliplr(t)],[curve1'; fliplr(curve2)'], ...
        C(i_set,:),'FaceAlpha',.3,'EdgeColor','none');
    hp.Annotation.LegendInformation.IconDisplayStyle = 'off';
    plot(t,ymean,'Color',C(i_set,:),'LineWidth',2,'DisplayName',['training set ',num2str(i_set)]);
end
u = legendUnq();
legend(u,'Box','off','color','none','Location','best','EdgeColor','w','TextColor','w')
set(gca, 'color', 'none', 'XColor','w', 'YColor','w', 'ZColor','w');
grid off
axis tight
set(gcf, 'color', 'none'); 
set(gcf, 'Position', [50 50 700 200]);
xticks(1:ntrials); xticklabels(stims(X));
ylabel('% hits')

pathout = fullfile('W:\scratch\gfriedri\caudtomm\ev_data\experiments\figures\');
FileOut = fullfile(pathout,strcat(classifier,'classified_by_trainset','.fig'));
if s; savefig(h,FileOut); end
FileOut = fullfile(pathout,strcat(classifier,'classified_by_trainset','.svg'));
if s; plot2svg(FileOut,h); end


novel_stims = {'Trp','Ser','Leu'};
familiar_stims = {'Arg','Ala','His'};
mask_novel = double(ismember(TestLabs_MAT, novel_stims));
mask_novel = permute(repmat(mask_novel,[1,1,nfish]),[3,1,2]);
mask_novel(mask_novel==0) = nan;
mask_familiar = double(ismember(TestLabs_MAT, familiar_stims));
mask_familiar = permute(repmat(mask_familiar,[1,1,nfish]),[3,1,2]);
mask_familiar(mask_familiar==0) = nan;


totFractionCorrectLabAllStims = squeeze(nanmean(correctTest,2));
totFractionCorrectLabSH = squeeze(nanmean(correctTestSH,2));
totFractionCorrectLabNovelStims = squeeze(nanmean(correctTest.*mask_novel,2));
totFractionCorrectLabFamiliarStims = squeeze(nanmean(correctTest.*mask_familiar,2));

totFractionCorrectLab = totFractionCorrectLabAllStims;
offset = totFractionCorrectLab(:,1);
totFractionCorrectLab = totFractionCorrectLab-offset;
maxval = mean(totFractionCorrectLab(:,3:5),2);
totFractionCorrectLab = totFractionCorrectLab./maxval;


%%
lip = .4;
h= figure; hold on
ylim([-.1,1.1]), xlim([1-lip, nsets+lip])
csh = [.6 .6 .6];
chancelv = 1/nstims;
line(xlim,repelem(chancelv,2), ...
    'Color','c','LineWidth',3,'LineStyle',':','DisplayName','chance')
for i_set = 1:nsets
    x = i_set+[-1, 1]*lip*2/3;
    
    curve1 = repelem(max(ylim),2); curve2 = repelem(min(ylim),2);
    hp = patch([x,fliplr(x)],[curve1'; fliplr(curve2)'], ...
        C(i_set,:),'FaceAlpha',.3,'EdgeColor','none');
    hp.Annotation.LegendInformation.IconDisplayStyle = 'off';

    line(x,repelem(nanmean(totFractionCorrectLabSH(:,i_set)),2), ...
        'Color',csh,'LineWidth',5,'DisplayName','shuffle')
end

y = totFractionCorrectLab;
jit = .1*randn(size(y));
x = repmat(1:nsets,[size(y,1),1]) + jit;
scatter(x, y,30,'w','filled')
plot(x',y','w','DisplayName','data')




m = [];
for i_set = 1:nsets
    x = i_set+[-1, 1]*lip*2/3;

    y = nanmean(totFractionCorrectLab(:,i_set));
    line(x,repelem(y,2),'Color','r','LineWidth',5)

    m = [m,y];
end

err = std(totFractionCorrectLab);
errorbar(m,err,'r');

u = legendUnq();
% legend(u,'Box','off','color','none','Location','best','EdgeColor','w','TextColor','w')

xticks(1:nsets)
xlabel('template trial #')
ylabel('overall performance on test data')
grid off
set(gca, 'color', 'none', 'XColor','w', 'YColor','w', 'ZColor','w');
set(gcf, 'color', 'none'); 
set(gcf, 'Position', [50 50 200 400])
