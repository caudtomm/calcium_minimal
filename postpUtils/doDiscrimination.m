function out = doDiscrimination(data, labs, varargin)
% doDiscrimination performs block-based multiclass decoding of neural
% population activity using a template-matching classifier, with optional
% within-trial neuron shuffling to estimate a chance-level baseline.
%
% INPUTS:
%   data  - numeric array of neural activity with one of the following shapes:
%             [cells x trials]
%             [time x cells x trials]
%           If time is present, activity is averaged across the time dimension
%           using mean(...,'omitmissing').
%
%   labs  - cell array [trials x 1] of trial labels. Label identity and
%           repetition structure are inferred from order of appearance.
%
% NAME-VALUE PAIRS (optional):
%   'method'         - distance metric passed to pdist for template matching
%                      (default: 'correlation')
%
%   'trainblockmode' - training block definition based on label repetitions:
%                        'single'  : train on one repetition index
%                        '3blocks' : train on sliding windows of three
%                                    consecutive repetition indices
%                      (default: 'single')
%
%   'classifier'     - decoding method. Currently supported:
%                        'template_match' (default)
%                      Other options (e.g. 'SVM') are present but not functional.
%
%   'nshuffles'      - number of shuffle iterations for baseline estimation
%                      (default: 50)
%
% CLASSIFICATION PROCEDURE:
%   - Trials are grouped into repetition blocks separately for each label.
%   - Training trials are selected based on repetition index; all remaining
%     trials are used for testing.
%   - For each class, a template is computed as the mean activity vector
%     across its training trials.
%   - Each trial is assigned the label of the nearest template according to
%     the chosen distance metric.
%
% SHUFFLE BASELINE:
%   - For each shuffle iteration, neuron identities are randomly permuted
%     independently within each trial, preserving trial-wise activity
%     distributions.
%
% OUTPUT:
%   out - structure containing decoding results:
%       .input_labs              - original trial labels
%       .train_trials            - [nTrials x nSets] logical training mask
%       .test_trials             - [nTrials x nSets] logical test mask
%       .predicted_labs          - [nTrials x nSets] predicted labels
%       .predicted_labs_SH       - [nTrials x nSets x nShuffles] shuffled predictions
%       .prediction_iscorrect    - correctness matrix (1/0/NaN)
%       .prediction_iscorrect_SH - shuffled correctness matrix
%       .prediction_confidence   - relative distance margin between best and
%                                  second-best template (not a probability)
% Notes:
% - Trials are grouped into repetition blocks based on label repetition order.
% - Shuffle testing permutes the neuron identity independently for each trial.
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

repetitions = 1:max(label_repetitions);
switch trainblockmode % # TODO : this should be adaptive to the num of repetitions of each label found in the data
    case '3blocks'
        trainblocksets = arrayfun(@(x) x:x+2, 1:repetitions(end)-2, 'UniformOutput', false);
    case 'single'
        trainblocksets = num2cell(repetitions);
    otherwise
        error('specified training blocks mode is unknown')
end

nsets = numel(trainblocksets);

% initialize output
out.input_labs = labs;
out.train_trials = false(ntrials,nsets);
out.test_trials = false(ntrials,nsets);
out.predicted_labs = cell(ntrials,nsets);
out.predicted_labs_SH = cell(ntrials,nsets,nshuffles);
out.prediction_iscorrect = nan(ntrials,nsets);
out.prediction_iscorrect_SH = nan(ntrials,nsets,nshuffles);
out.prediction_confidence = nan(ntrials,nsets);


%% action
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

        switch lower(classifier)
            case "svm"
                res = fit_SVM(trainData',trainlabs,testData',stims);
            case "template_match"
                res = template_matching(trainData', trainlabs, testData', stims, method);
            case 'lda'
                res = lindiscrim(trainData',trainlabs,testData',stims,'linear'); % needs multiple train-examples per class
            case 'qda'
                res = lindiscrim(trainData',trainlabs,testData',stims,'quadratic'); % needs multiple train-examples per class
            case 'dbd'
                res = dbd(trainData',trainlabs,testData',stims);
            otherwise
                error('unknown classifier')
        end

        % combine predictions for training and testing trials into one cell
        % array (don't worry, we are holding onto the indices)
        predicted_labs = cell(ntrials,1);
        predicted_labs(trials_train) = res.predictions_trainData;
        predicted_labs(trials_test) = res.predictions_testData;

        % combine confidence values
        prediction_confidence = nan(ntrials,1);
        prediction_confidence(trials_train) = res.confidence_trainData;
        prediction_confidence(trials_test) = res.confidence_testData;

        % make an array to note only whether the predictions were correct
        prediction_iscorrect = double(cellfun(@isequal, labs, predicted_labs));
        missing_data = cellfun(@isempty, predicted_labs);
        prediction_iscorrect(missing_data)=nan;

        switch i
            case 1 % non-shuffled data
                out.train_trials(:,i_set) = trials_train;
                out.test_trials(:,i_set) = trials_test;
                out.predicted_labs(:,i_set) = predicted_labs;
                out.prediction_iscorrect(:,i_set) = prediction_iscorrect;
                out.prediction_confidence(:,i_set) = prediction_confidence;

            otherwise % shuffled data
                out.predicted_labs_SH(:,i_set,i-1) = predicted_labs;
                out.prediction_iscorrect_SH(:,i_set,i-1) = prediction_iscorrect;
        end
        
        
    end
end


end


%% Functions

% function [yfit, predictions] = fit_SVM(trainData,trainlabs,testData,stims)
%     [svm, accuracy, predictions] = trainSVM(trainData,trainlabs,stims);
%     yfit = svm.predictFcn(testData);
% end

function out = fit_SVM(trainData, trainlabs, testData, stims)
    % ensure categorical labels for ECOC
    trainlabs_cat = categorical(trainlabs, stims);

    svm = fitcecoc(trainData, trainlabs_cat, ...
        'Coding','onevsall', ...
        'Learners','linear', ...
        'Verbose',0);

    [pred_train, score_train] = predict(svm, trainData);
    [pred_test,  score_test ] = predict(svm, testData);

    out.predictions_trainData = cellstr(pred_train);
    out.predictions_testData  = cellstr(pred_test);

    out.confidence_trainData = scoreMargin(score_train);
    out.confidence_testData  = scoreMargin(score_test);

    function conf = scoreMargin(scores)
        % scores: [nSamples x nClasses]
        scores = sort(scores,2,'descend');
        conf = (scores(:,1) - scores(:,2)) ./ abs(scores(:,1));
    end
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

function out = train_RNN(trainData, trainlabs, testData, stims)
    % Placeholder for future RNN implementation
end

function out = lindiscrim(trainData, trainlabs, testData, stims, discrimType)
    trainlabs_cat = categorical(trainlabs, stims);
    % Ensure there are multiple entries for each category in trainlabs (required for this classifier)
    if any(histcounts(trainlabs_cat) < 2)
        error('Each category must have at least two training set entries for lda/qda.');
    end

    switch discrimType
        case 'linear'
            gamma = 0;
        case 'quadratic'
            gamma = 1;
        otherwise
            error('unknown discrimType for lindiscrim')
    end

    lda = fitcdiscr(trainData, trainlabs_cat, ...
        'DiscrimType',discrimType, ... % lda if 'linear', qda is 'quadratic'
        'Gamma', gamma); 

    [pred_train, score_train] = predict(lda, trainData);
    [pred_test,  score_test ] = predict(lda, testData);

    out.predictions_trainData = cellstr(pred_train);
    out.predictions_testData  = cellstr(pred_test);

    out.confidence_trainData = scoreMargin(score_train);
    out.confidence_testData  = scoreMargin(score_test);

    function conf = scoreMargin(scores)
        scores = sort(scores,2,'descend');
        conf = (scores(:,1) - scores(:,2)) ./ abs(scores(:,1));
    end
end

function out = dbd(trainData, trainlabs, testData, stims)
    nstims = numel(stims);
    nvars  = size(trainData,2);

    bases = nan(nstims, nvars);
    for i = 1:nstims
        bases(i,:) = mean(trainData(ismember(trainlabs,stims{i}),:),1,'omitmissing');
    end

    [out.predictions_trainData, out.confidence_trainData] = ...
        predictLabels(trainData, bases, stims);

    [out.predictions_testData, out.confidence_testData] = ...
        predictLabels(testData, bases, stims);

    function [predictions, confidence] = predictLabels(data, bases, labels)
        scores = data * bases'; % projection
        [scores, idx] = sort(scores,2,'descend');

        predictions = labels(idx(:,1));
        confidence  = (scores(:,1) - scores(:,2)) ./ abs(scores(:,1));
    end
end
