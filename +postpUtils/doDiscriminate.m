function [totFractionCorrectLab,totFractionCorrectLabSH, m] = doDiscriminate(experiment, classifier, stims2use, trainblockmode)

todo_fish = find(ismember(experiment.summaryTable.group,{'trained1';'trained2';'trained1 (T-R-S-H-A-ACSF/L)';'uncoupled'}));% options : {'previousnaive';'naïve';'trained1';'trained2';'trained1 (T-R-S-H-A-ACSF/L)';'uncoupled'};

nfish = numel(todo_fish);
if isempty(stims2use); stims2use = {'Trp','Ser','Ala','Food'}; end
nstims = numel(stims2use);
trialn2usein = 1:5;
if isempty(s); s = 0; end
method = 'cosine';

ntrials=30;

oneblock = [1:nstims];
switch trainblockmode
    case '3blocks'
        trainblocksets = {1:3, 2:4, 3:5};
    case 'single'
        trainblocksets = {1,2,3,4,5};
    otherwise
        error('specified training blocks mode is unknown')
end

nsets = numel(trainblocksets);

nshuffles = 50;

if isempty(classifier); classifier = "template_match"; end

% time period to include (relative to [stim_on, stim_on) [s]
if ~exist('tframe','var') || isempty(tframe); tframe = [1 20]; end
tframe = [1 20]

%% action
correctLab = [];
correctTest = [];
correctLabSH = [];
correctTestSH = [];
TestLabs_MAT = [];
for i_fish = 1:nfish
    % needs: data from multiple fish, 
    [a, X, stims] = selectData(experiment, todo_fish(i_fish),stims2use,trialn2usein);
    
    ncells = size(a,2);
    ntrials = numel(X);
    L = size(a,1)/ntrials;
    
    a = traceFormat(a,L);
    a = squeeze(nanmean(a,1));
    
    for i_set = 1:nsets
        
        trials_train = [];
        for i_stim = 1:nstims
            thisstim = find(ismember(stims,stims2use(i_stim)));
            thisstimidx = find(ismember(X,thisstim));
            if isempty(thisstimidx); continue;end
            trials_train = [trials_train; thisstimidx(trainblocksets{i_set})];
        end

        trials_test = [1:ntrials]'; trials_test(trials_train) = [];
    
        trainlabs = stims(X(trials_train))';
        testlabs = stims(X(trials_test))';

        if i_fish == 1; TestLabs_MAT = [TestLabs_MAT, testlabs]; end

        for i = 1:nshuffles+1 % 1x data + 50x shuffle
            switch i
                case 1
                    tmp = a;
                otherwise
                    a_shuf = a;
                    for i_trial = 1:ntrials
                        a_shuf(i_trial,:) = a_shuf(i_trial,randperm(ncells));
                    end
                    tmp = a_shuf;
            end
%             tmp = tmp - nanmean(tmp,1); % men-subtract each feature (neuron)
            trainData = tmp(trials_train,:);
            testData = tmp(trials_test,:);
    
            switch classifier
                case "SVM"
                    yfit = fit_SVM(trainData,trainlabs,testData,stims2use);
                case "template_match"
                    yfit = template_matching(trainData, trainlabs, testData, stims2use, method);
                otherwise
                    error('unknown classifier')
            end
    
            correct_test = cellfun(@isequal, testlabs, yfit);
            correct_train = cellfun(@isequal, trainlabs, predictions);
            thiscorrect = nan(ntrials,1);
            thiscorrect(trials_test) = correct_test;
            thiscorrect(trials_train) = correct_train;

            switch i
                case 1
                    try
                        correctLab = [correctLab, thiscorrect];
                        correctTest = [correctTest, correct_test];
                    catch
                        thiscorrect = [thiscorrect; nan(size(correctLab,1)-length(thiscorrect),1)];
                        correct_test = [correct_test; nan(size(correctTest,1)-length(correct_test),1)];

                        correctLab = [correctLab, thiscorrect];
                        correctTest = [correctTest, correct_test];
                    end
                otherwise
                    try
                        correctLabSH = [correctLabSH, thiscorrect];
                        correctTestSH = [correctTestSH, correct_test];

                    catch
                        thiscorrect = [thiscorrect; nan(size(correctLabSH,1)-length(thiscorrect),1)];
                        correct_test = [correct_test; nan(size(correctTest,1)-length(correct_test),1)];
    
                        correctLabSH = [correctLabSH, thiscorrect];
                        correctTestSH = [correctTestSH, correct_test];
                    end
            end
            
            
        end
    end
end

correctLab = permute(traceFormat(correctLab',nsets),[2,3,1]);
correctTest = permute(traceFormat(correctTest',nsets),[2,3,1]);

correctLabSHtmp = permute(traceFormat(correctLabSH',nsets*nshuffles),[2,3,1]);
correctLabSH = nan(size(correctLab));
for i_set = 1:nsets
    idx = (i_set-1)*nshuffles + 1 : i_set*nshuffles;
    correctLabSH(:,:,i_set) = nanmean(correctLabSHtmp(:,:,idx),3);
end
clear correctLabSHtmp

correctTestSHtmp = permute(traceFormat(correctTestSH',nsets*nshuffles),[2,3,1]);
correctTestSH = nan(size(correctTest));
for i_set = 1:nsets
    idx = (i_set-1)*nshuffles + 1 : i_set*nshuffles;
    correctTestSH(:,:,i_set) = nanmean(correctTestSHtmp(:,:,idx),3);
end
clear correctTestSHtmp

% return
out.correctLab = correctLab;
out.correctTest = correctTest;
out.correctLabSH = correctLabSH;
out.correctTestSH = correctTestSH;

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

FileOut = fullfile(pathout,strcat(classifier,'classified_by_trainset2','.fig'));
if s; savefig(h,FileOut); end
FileOut = fullfile(pathout,strcat(classifier,'classified_by_trainset2','.svg'));
if s; plot2svg(FileOut,h); end

[p,h,stats] = ranksum(totFractionCorrectLab(:,1),totFractionCorrectLab(:,end))
y = totFractionCorrectLab(:,end-2:end);
[p,h,stats] = ranksum(totFractionCorrectLab(:,1),y(:))

end


%% Functions

function yfit = fit_SVM(trainData,trainlabs,testData,stims2use)
    [svm, accuracy, predictions] = pospUtils.trainSVM(trainData,trainlabs,stims2use);
    yfit = svm.predictFcn(testData);
end

function yfit = template_matching(trainData, trainlabs, testData, stims2use, method)
    [templates, templates_weights] = deal(nan(nstims, size(trainData,2)));
    for i_stim = 1:nstims
        thisstim_trials = ismember(trainlabs,stims2use{i_stim});
        templates(i_stim,:) = nanmean(trainData(thisstim_trials,:),1);
%                         tmp = std(trainData(thisstim_trials,:),1);
%                         tmp = nanmin(tmp) ./ tmp;
%                         templates_weights(i_stim,:) = tmp;
    end

    %method = 'cosine';
    distances = nan(numel(trials_train),nstims);
    for i_trial = 1:numel(trials_train)
        thistrial = trainData(i_trial,:);
        for i_template = 1:nstims
            distances(i_trial,i_template) = ...
                pdist([templates(i_template,:); thistrial],method);
        end
    end
    [distances, idx] = sort(distances,2);
    predictions = stims2use(idx(:,1))';
    predictions_confidence = diff(distances(:,1:2),[],2) ./ 2; % DIVISION BY 2 IN THE CASE OF COSINE DISTANCE


    distances = nan(numel(trials_test),nstims);
    for i_trial = 1:numel(trials_test)
        thistrial = testData(i_trial,:);
        for i_template = 1:nstims
            distances(i_trial,i_template) = ...
                pdist([templates(i_template,:); thistrial],method);
        end
    end
    [distances, idx] = sort(distances,2);
    yfit = stims2use(idx(:,1))';
    yfit_confidence = diff(distances(:,1:2),[],2);
end
