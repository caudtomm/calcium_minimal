function out = doDiscriminate(data, labs, varargin)
arguments
    data double % [variables x samples] or [cells x trials]; or [time x variables x samples]
    labs cell % {samples x 1} or {trials x 1} of strings
end
arguments (Repeating)
    varargin
end

% Set default values
method = 'correlation';
trainblockmode = 'single';
classifier = 'template_match'; % 'SVM' or 'template_match'
nshuffles = 50;

% Parse name-value pairs
if ~isempty(varargin)
    for k = 1:2:length(varargin)
        switch lower(varargin{k})
            case 'method'
                method = varargin{k+1};
            case 'trainblockmode'
                trainblockmode = varargin{k+1};
            case 'classifier'
                classifier = varargin{k+1};
            case 'nshuffles'
                nshuffles = varargin{k+1};
        end
    end
end

if ndims(data)==3 % [time x variables x samples]
    % compress time to average
    data = squeeze(mean(data,'omitmissing'));
end
[ncells, ntrials] = size(data);
stims = unique(labs);
nstims = numel(stims);

% # TODO : validation of input data and labels (same number of trials, same number of repetitions across labels, etc.)

% Compute repetition number for each label (e.g. A,B,C,A,B,C -> 1,1,1,2,2,2)
label_repetitions = zeros(size(labs));
for i = 1:nstims
    idx = find(ismember(labs, stims{i}));
    label_repetitions(idx) = 1:numel(idx);
end

repetitions = 1:5;
oneblock = [1:nstims];
switch trainblockmode % # TODO : this should be adaptive to the num of repetitions of each label found in the data
    case '3blocks'
        trainblocksets = {1:3, 2:4, 3:5};
    case 'single'
        trainblocksets = {1,2,3,4,5};
    otherwise
        error('specified training blocks mode is unknown')
end

nsets = numel(trainblocksets);

%% action
correctLab = [];
correctTest = [];
correctLabSH = [];
correctTestSH = [];


for i_set = 1:nsets
    % specify trial indices to train and test on
    trials_train = ismember(label_repetitions, trainblocksets{i_set});
    trials_test = ~trials_train;

    % correct labels for training and testing
    trainlabs = labs(trials_train);
    testlabs = labs(trials_test);

    for i = 1:nshuffles+1 % 1x data + 50x shuffle
        switch i
            case 1
                tmp = data;
            otherwise
                a_shuf = data;
                for i_trial = 1:ntrials % for loop to shuffle each trial independently
                    a_shuf(:,i_trial) = a_shuf(randperm(ncells),i_trial);
                end
                tmp = a_shuf;
        end
        trainData = tmp(:,trials_train);
        testData = tmp(:,trials_test);

        switch classifier
            case "SVM"
                % [yfit, predictions] = fit_SVM(trainData',trainlabs,testData',stims);
            case "template_match"
                out = template_matching(trainData', trainlabs, testData', stims, method);
            otherwise
                error('unknown classifier')
        end

        correct_test = double(cellfun(@isequal, testlabs, out.predictions_testData));
        missing_data = cellfun(@isempty, out.predictions_testData);
        correct_test(missing_data)=nan;

        correct_train = double(cellfun(@isequal, trainlabs, out.predictions_trainData));
        missing_data = cellfun(@isempty, out.predictions_trainData);
        correct_train(missing_data)=nan;

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
            otherwise % # TODO this should make a 3d matrix of [trials x sets x nshuffles]
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


%[p,h,stats] = ranksum(totFractionCorrectLab(:,1),totFractionCorrectLab(:,end))
% y = totFractionCorrectLab(:,end-2:end);
% [p,h,stats] = ranksum(totFractionCorrectLab(:,1),y(:))

end


%% Functions

function [yfit, predictions] = fit_SVM(trainData,trainlabs,testData,stims)
    [svm, accuracy, predictions] = postpUtils.trainSVM(trainData,trainlabs,stims);
    yfit = svm.predictFcn(testData);
end

function out = template_matching(trainData, trainlabs, testData, stims, method)
    nstims = numel(stims);
    ncells = size(testData,2);

    % define templates by averaging across 'training' vectors with the same label
    % templates follow the same order as stims
    templates = nan(nstims, ncells);
    for i_stim = 1:nstims
        thisstim_trials = ismember(trainlabs,stims{i_stim});
        templates(i_stim,:) = mean(trainData(thisstim_trials,:),1,'omitmissing');
    end

    [predictions_trainData, confidence_trainData] = predictLabels(trainData,templates,stims);
    [predictions_testData, confidence_testData] = predictLabels(testData,templates,stims);

    out.predictions_trainData = predictions_trainData;
    out.predictions_testData = predictions_testData;
    out.confidence_trainData = confidence_trainData;
    out.confidence_testData = confidence_testData;

    function [predictions, confidence] = predictLabels(data,templates,template_labels)
        [n_samples,n_variables] = size(data);
        n_templates = size(templates,1);
        distances = nan(n_samples,n_templates);
        for i_sample = 1:n_samples
            thissample = data(i_sample,:);
            for i_template = 1:n_templates
                distances(i_sample,i_template) = ...
                    pdist([templates(i_template,:); thissample],method);
            end
        end
        [distances, idx] = sort(distances,2);
        existing_data = ~isnan(distances(:,1));
    
        predictions = cell(n_samples,1);
        confidence = nan(n_samples,1);
        predictions(existing_data) = template_labels(idx(existing_data,1));
        confidence(existing_data) = ...
            (distances(existing_data,2) - distances(existing_data,1)) ./ (distances(existing_data,2));
    end
end
