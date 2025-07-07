function out = plotDiscrimination(predictions,plotType,cfg,varargin)
% plotDiscrimination visualizes classification results produced by doDiscriminate.
% This low-level function is reusable for displaying performance metrics such as
% accuracy, confidence, and shuffle comparisons across different datasets or subjects.
% Assumes equal representation for each provided label (stimulus)
%
% INPUTS:
%   predictions - cell array [nSubjects x 1], each cell containing a struct output
%                 from doDiscriminate. Each struct must include the following fields:
%
%       .input_labs               - {nTrials x 1} original labels for each trial
%       .train_trials             - [nTrials x nSets] logical array indicating training trials
%       .test_trials              - [nTrials x nSets] logical array indicating test trials
%       .predicted_labs           - {nTrials x nSets} predicted labels for original data
%       .predicted_labs_SH        - {nTrials x nSets x nShuffles} predicted labels for shuffled data
%       .prediction_iscorrect     - [nTrials x nSets] array (0/1/NaN) of prediction correctness
%       .prediction_iscorrect_SH  - [nTrials x nSets x nShuffles] array of correctness for shuffled data
%       .prediction_confidence    - [nTrials x nSets] prediction confidence scores
%
%   plotType    - string specifying what to visualize (default: 'full').
%                 Supported values (implementation-dependent):
%                   'full'         - plot all available metrics
%                   'accuracy'     - only classification accuracy (real vs shuffled)
%                   'confidence'   - only confidence values
%                   'perStimulus'  - per-label accuracy breakdown
%
%   cfg         - PlotConfig object with user-defined visual settings
%
% OUTPUT:
%   out - (optional)
%
% USAGE:
%   plotDiscrimination({res1; res2}, 'accuracy', myCfg);
%
% NOTES:
% - This function is designed for modular use in larger pipelines (e.g., group analysis, reports).
% - Predictions from different subjects/datasets must follow the same format
%   (see output of postpUtils.doDiscriminate).

arguments
    predictions cell % cell array [nsubjects 1] of double [t, cells, trials] (sorted!)
    plotType char = 'full'
    cfg = PlotConfig()
end
arguments (Repeating)
    varargin
end

%
ntrials = cellfun(@(x) numel(x.input_labs),predictions);
% Set default values
focus_trials = arrayfun(@(n) (1:n)', ntrials, 'UniformOutput', false);
do_zscore = false;

% Parse name-value pairs
if ~isempty(varargin)
    for k = 1:2:length(varargin)
        switch lower(varargin{k})
            case 'focustrials'
                focus_trials = varargin{k+1};
            case 'zscore'
                do_zscore = varargin{k+1};
        end
    end
end

% from input predictions, only keep focus trials
nsubjects = numel(predictions);
p = cell(nsubjects,1);
for i=1:nsubjects
    p{i}.input_labs = predictions{i}.input_labs(focus_trials{i});
    p{i}.train_trials = predictions{i}.train_trials(focus_trials{i},:);
    p{i}.test_trials = predictions{i}.test_trials(focus_trials{i},:);
    p{i}.predicted_labs = predictions{i}.predicted_labs(focus_trials{i},:);
    p{i}.predicted_labs_SH = predictions{i}.predicted_labs_SH(focus_trials{i},:,:);
    p{i}.prediction_iscorrect = predictions{i}.prediction_iscorrect(focus_trials{i},:);
    p{i}.prediction_iscorrect_SH = predictions{i}.prediction_iscorrect_SH(focus_trials{i},:,:);
    p{i}.prediction_confidence = predictions{i}.prediction_confidence(focus_trials{i},:);
end

% useful metrics (based on first subject)
[ntrials,nsets,nshuffles] = size(p{1}.predicted_labs_SH);
stims = unique(p{1}.input_labs);
nstims = numel(stims);


%% Plot onto provided axes

switch plotType
    case 'performance_lines' % Plot all in one big matrix
        out = performanceLines();
    otherwise
        error('Requested plot type is unknown.')
end

%% plotting
function out = plot1()
    % for i_set = 1:nsets
    %     subplot(2,nsets,i_set); imagesc(correctLab(:,:,i_set))
    %     title(['training set ',num2str(i_set)])
    %     ylabel('fish')
    % end
    % for i_set = 1:nsets
    %     subplot(2,nsets,i_set+nsets); imagesc(correctLabSH(:,:,i_set))
    % end
    % xticks(1:numel(X)); xticklabels(stims(X));
end

function out = plot2()
    % C = optimalcolors(nsets+1); C = C(2:end,:);
    % 
    % h = figure; hold on
    % t = 1:ntrials;
    % for i_set = 1:nsets
    %     tmp = correctLab(:,:,i_set);
    %     ymean = 100* nanmean(tmp,1);
    %     ystd = 100* squeeze(std(tmp,[],1,'omitnan'));
    % 
    %     curve1 = ymean + ystd; curve2 = ymean - ystd;
    %     hp = patch([t,fliplr(t)],[curve1'; fliplr(curve2)'], ...
    %         C(i_set,:),'FaceAlpha',.3,'EdgeColor','none');
    %     hp.Annotation.LegendInformation.IconDisplayStyle = 'off';
    %     plot(t,ymean,'Color',C(i_set,:),'LineWidth',2,'DisplayName',['training set ',num2str(i_set)]);
    % end
    % u = legendUnq();
    % legend(u,'Box','off','color','none','Location','best','EdgeColor','w','TextColor','w')
    % set(gca, 'color', 'none', 'XColor','w', 'YColor','w', 'ZColor','w');
    % grid off
    % axis tight
    % set(gcf, 'color', 'none'); 
    % set(gcf, 'Position', [50 50 700 200]);
    % xticks(1:ntrials); xticklabels(stims(X));
    % ylabel('% hits')
end

function out = performanceLines()
    out = [];

    hold on

    if ~do_zscore
        lip = .4;
        ylim([-.1,1.1]), xlim([1-lip, nsets+lip])
        csh = [.6 .6 .6]; % color for shuffled data
    
        % chance level dotted line # BUG this assumes equal representation of all labels 
        chancelv = 1/nstims;
        line(xlim,repelem(chancelv,2), ...
            'Color','c','LineWidth',3,'LineStyle',':','DisplayName','chance')
    
        %% shuffled data performance averages
        % get sh performance by averaging over all the fish and trials for each
        % training set
        performanceSH = cell2mat(cellfun(@(x) x.prediction_iscorrect_SH, p, ...
            'UniformOutput',false)); % [all_trials x nsets x nshuffles]
        performanceSH = mean(performanceSH,[1,3],'omitmissing');
        % plot average performances as little line segments
        x = [1:nsets]'+[-1, 1]*lip*2/3;
        line(x',repmat(performanceSH,2,1), ...
                'Color',csh,'LineWidth',5,'DisplayName','shuffle')
    end
    %% single-subject lines
    performance = cell2mat(cellfun(@(x) mean(x.prediction_iscorrect,"omitmissing"), p, ...
        'UniformOutput',false)); % [subjects x nsets]
    if do_zscore; performance = nanzscore(performance,[],2); end
    jit = .1*randn(size(performance));
    x = repmat(1:nsets,[nsubjects,1]) + jit;
    scatter(x, performance,30,cfg.axcol,'filled')
    plot(x',performance','Color',cfg.axcol,'DisplayName','data')

    %% average performance over subjects
    x = 1:nsets;
    y = mean(performance,'omitmissing');
    line(x,y,'Color','r','LineWidth',5)
    err = std(performance);
    errorbar(y,err,'r');

    %% cosmetics
    u = legendUnq();
    legend(u,'Box','off','color','none','Location','best','EdgeColor',cfg.axcol,'TextColor',cfg.textcol)

    xticks(1:nsets)
    xlabel('template trial #')
    ylabel('overall performance on test data')
    grid off
    set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
    set(gcf, 'color', cfg.bgcol); 
    % set(gcf, 'Position', [50 50 200 400])

    %% return
    out = performance;
end

end