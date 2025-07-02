function out = doDiscriminate(data, labs, varargin)
% doDiscriminate performs classification on neural response data using a specified
% classifier and cross-validation scheme, and optionally computes a shuffle baseline.
%
% INPUTS:
%   data      - double matrix of neural activity, with one of the following shapes:
%                 [variables x trials]        (e.g., cells x trials)
%                 [time x variables x trials] (will be averaged across time)
%   labs      - cell array {trials x 1} of string labels, one per trial
%
% NAME-VALUE PAIRS (optional):
%   'method'          - similarity metric for classification (default: 'correlation')
%   'trainblockmode'  - cross-validation scheme: 'single' or '3blocks' (default: 'single')
%   'classifier'      - classification method: 'template_match' (default) or 'SVM' (not implemented)
%   'nshuffles'       - number of shuffle iterations for significance testing (default: 50)
%
% OUTPUT:
%   out - structure containing classification results:
%       .input_labs               - original trial labels
%       .train_trials             - [nTrials x nSets] logical matrix of training indices
%       .test_trials              - [nTrials x nSets] logical matrix of testing indices
%       .predicted_labs           - [nTrials x nSets] cell array of predicted labels
%       .predicted_labs_SH        - [nTrials x nSets x nShuffles] cell array of shuffled predictions
%       .prediction_iscorrect     - [nTrials x nSets] array of correctness (1=correct, 0=incorrect, NaN=missing)
%       .prediction_iscorrect_SH  - [nTrials x nSets x nShuffles] array of shuffled correctness
%       .prediction_confidence    - [nTrials x nSets] array of prediction confidence values
%
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

        switch classifier
            % case "SVM"
                % [yfit, predictions] = fit_SVM(trainData',trainlabs,testData',stims);
            case "template_match"
                res = template_matching(trainData', trainlabs, testData', stims, method);
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
