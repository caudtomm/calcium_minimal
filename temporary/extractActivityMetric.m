function [output_struct, datatab, datatab_baseline] = extractActivityMetric(experiment,params)
disp(' ')

% initialize variables
time_range = params.time_range;
time_range_baseline = [-diff(params.time_range),0]-2; % [from stim_on, from stim_on] (sec)
learned_stim_tag = {'all CS+';
                    'all CS-';
                    'all familiar';
                    'all novel'};
switch params.dimension
    case 'over_time' % yields one value per cell, per trial
        desired_dimension = 1;
    case 'over_cells' % yields one value per frame, per trial
        desired_dimension = 2;
    otherwise
        error('desired dimension is invalid.')
end

% select which fish are part of this analysis
todo_fish = select_fish(experiment, params.groups);

% initialize loop outputs
fish_stats = [];
outputMAT = []; outputMAT_baseline = [];


% for each fish, extract a single vector per trial
for i_fish = 1:numel(todo_fish)
    thisfish = todo_fish(i_fish);
    fprintf('fish num : %s\n', num2str(thisfish))
    data = experiment.series{thisfish}.data;

    % hacky shit
    if iscolumn(data.stim_type)
        data.stim_type = data.stim_type';
    end

    output_vector = []; output_vector_baseline = [];
    
    % extract fish-specific vars
    fs = data.meta.framerate;
    
    % define desired time-interval
    interval = floor( (time_range(1) + data.stim_on_sec) * fs ) : ...
               floor( (time_range(2) + data.stim_on_sec) * fs );
    interval_baseline = floor( (time_range_baseline(1) + data.stim_on_sec) * fs ) : ...
                        floor( (time_range_baseline(2) + data.stim_on_sec) * fs );
    fprintf('interval length check : %s\n', string(length(interval_baseline)==length(interval)))

    % define desired trials
    if ismember(params.stimuli,learned_stim_tag) % dynamic stimulus group (from learning)
        desired_stim_type = select_trained_stimuli(experiment.series{thisfish}.group,params.stimuli);
    else
        desired_stim_type = params.stimuli;
    end
    relative_trial_num = get_relative_trial_num(data.stim_type); 
    desired_trials = ismember(data.stim_type,desired_stim_type) & ismember(relative_trial_num,params.todo_reltrialnum);
    stim_type = data.stim_type(desired_trials);
    trial_num = data.trial_num(desired_trials);
    relative_trial_num = relative_trial_num(desired_trials);

    % extract activity traces
    tracestmp = traceFormat(data.traces,data.L); % time x cells x trials

    % define desired cells
    cells = 1:data.N;
    if params.only_top_variant_units
        FindTopUnits
        cells = alltopunits';
    end
    
    % select only needed traces
    traces = tracestmp(interval,cells,desired_trials);
    traces_baseline = tracestmp(interval_baseline,cells,desired_trials);

    if isempty(traces);disp('traces empty.');continue;end

    % get desired activity vector
    switch params.metric
        case 'max_activity'
            output_vector = squeeze(nanmax(traces,[],desired_dimension)); % (cells or time) x trials
            output_vector_baseline = squeeze(nanmax(traces_baseline,[],desired_dimension));
        case 'avg_activity'
            output_vector = squeeze(nanmean(traces,desired_dimension)); %  (cells or time) x trials
            output_vector_baseline = squeeze(nanmean(traces_baseline,desired_dimension));
        case 'variance'
            output_vector = squeeze(std(traces,[],desired_dimension,'omitnan')).^2; % (cells or time) x trials
            output_vector_baseline = squeeze(std(traces_baseline,[],desired_dimension,'omitnan')).^2;
        case 'lifetime_kurtosis'
            output_vector = squeeze(calculateSparseness(traces,'kurtosis-singletrials'));
            output_vector_baseline = squeeze(calculateSparseness(traces_baseline,'kurtosis-singletrials'));
        otherwise
            error('Requested metric invalid.')
    end

    % % normalize by baseline average
    % output_vector = output_vector ./ nanmean(output_vector_baseline,"all");
    % output_vector_baseline = output_vector_baseline ./ nanmean(output_vector_baseline,"all");

    % store loop outputs
    fish_stats{i_fish}.fish_id = thisfish;
    fish_stats{i_fish}.stim_type = stim_type;
    fish_stats{i_fish}.trial_num = trial_num;
    fish_stats{i_fish}.relative_trial_num = relative_trial_num;
    outputMAT{i_fish} = output_vector;
    outputMAT_baseline{i_fish} = output_vector_baseline;
    
end
output_struct.fish_stats = fish_stats;
output_struct.metricMAT = outputMAT;
output_struct.metricMATbaseline = outputMAT_baseline;
datatab = turnintoDataTable(output_struct.metricMAT,output_struct.metricMATbaseline,output_struct.fish_stats,params);
datatab_baseline = turnintoDataTable(output_struct.metricMATbaseline,output_struct.metricMATbaseline,output_struct.fish_stats,params);


end

%% Accessory functions

% Function to select fish to analyse, based on their groups
function todo_fish = select_fish(experiment, desired_groups)
    todo_fish = [];
    for i_fish = 1:numel(experiment.series)
        if ismember(experiment.series{i_fish}.group, desired_groups)
            todo_fish = [todo_fish;i_fish];
        end
    end
end

% # TODO : remove from here. this should only be in PlotConfig
% Function to extract group-specific identities of stimuli based on
% training status
function desired_stimuli = select_trained_stimuli(group,tag)
    all_groups = {'previousnaive';'naive';'trained1';'trained2';'trained1alt';'uncoupled'};
    all_CSplus = {{''};{''};{'Arg','Ala'};{'Ala','His'};{'Arg','Ala'};{''}};
    all_CSminus = {{''};{''};{'His'};{'Arg'};{'His'};{'Arg','Ala','His'}};
    all_familiar = {{''};{''};{'Arg','Ala','His'};{'Arg','Ala','His'};{'Arg','Ala','His'};{'Arg','Ala','His'}}; % CS+ and CS- together
    all_novel = {{'Trp','Ala','Ser','Food'};{'Trp','Ser','Leu'};{'Trp','Ser','Leu'};{'Trp','Ser','Leu'};{'Trp','Ser','Leu'};{'Trp','Ser','Leu'}};
    groups2stims = table(all_groups,all_CSplus,all_CSminus,all_familiar,all_novel);
    
    group_idx = ismember(groups2stims.all_groups,group);
    switch tag
        case 'all CS+'
            desired_stimuli = groups2stims.all_CSplus{group_idx};
        case 'all CS-'
            desired_stimuli = groups2stims.all_CSminus{group_idx};               
        case 'all familiar'
            desired_stimuli = groups2stims.all_familiar{group_idx};
        case 'all novel'
            desired_stimuli = groups2stims.all_novel{group_idx};
        otherwise
            error('stimulus group tag invalid.')
    end
end

% Function to turn stimulus identity train into relative trial numbers
function relative_trial_num = get_relative_trial_num(stim_type)
    relative_trial_num = [];
    if isempty(stim_type); return;end
    relative_trial_num = nan(1,length(stim_type));
    stims = unique(stim_type);
    for i_stim = 1:numel(stims)
        thisstim = stims{i_stim};
        thisstim_trials_idx = find(ismember(stim_type,thisstim));
        for i_trial = 1:numel(thisstim_trials_idx)
            relative_trial_num(thisstim_trials_idx(i_trial)) = i_trial;
        end
    end
end

%
function datatab = turnintoDataTable(dataMAT,baselineMAT,fish_stats,params)
    
    % set local variable names
    organizer = params.order_by;
    n_equals = params.n_equals;
    
    % initialize data table fields
    organizer_value = []; % field 1
    metric_value = []; % field 2
    
    % set data fields
    for i_fish = 1:numel(dataMAT)
    %     organizer_value = transpose(outputMAT{i_fish}.(organizer));
        thisdataMAT = dataMAT{i_fish}; % cells x trials
        thisbaselineMAT = baselineMAT{i_fish}; % cells x trials
        ncells = size(thisdataMAT,1);
        metric_value_fish = [];
        organizer_value_fish = [];
        
        % select cells according to their relative value in the chosen metric, if applicable
        if strcmp(params.dimension,'over_time')
            thisdataMAT = extract_quantile(thisdataMAT,thisbaselineMAT,params);
        end
    
        for i_trial = 1:size(thisdataMAT,2)
            metric_value_fish = [metric_value_fish; ...
                                thisdataMAT(:,i_trial)];
            organizer_value_fish = [organizer_value_fish; ...
                                   repelem(fish_stats{i_fish}.(organizer)(i_trial),ncells,1)];
        end
    
        switch n_equals
            case 'all'
            case 'fish'
                orgs = unique(organizer_value_fish);
    
                metric_tmp = nan(numel(orgs),1);
                for i_org = 1:numel(orgs)
                    idx = organizer_value_fish==orgs(i_org);
                    metric_tmp(i_org) = nanmean(metric_value_fish(idx),'all');
                end
    
                metric_value_fish = metric_tmp;
                organizer_value_fish = orgs;
            otherwise
                error('type of n is invalid.')
        end
        metric_value = [metric_value; metric_value_fish];
        organizer_value = [organizer_value; organizer_value_fish];
    end

    datatab = table(organizer_value,metric_value);
end


function dataMAT = extract_quantile(dataMAT,baselineMAT,params)
    % select indices
    switch params.isolate_quantile
        case 'all'
            idx = ones(size(dataMAT));
        case 'bottom'
            thisquantile = 1-params.quantile;
            quantileMSG(['bottom ',num2str(thisquantile)])
            idx = double(dataMAT<= quantile(baselineMAT,thisquantile,1));
        case 'top'
            thisquantile = params.quantile;
            quantileMSG(['top ',num2str(thisquantile)])
            idx = double(dataMAT>= quantile(baselineMAT,thisquantile,1));
        otherwise
            idx = ones(size(dataMAT));
            warning('Parameter property ''isolate_quantile'' lacks an associated function. All cells kept for further analysis');
    end
    idx(~idx)=nan;
    dataMAT = dataMAT.*idx;

    function quantileMSG(inset)
        msg = ['Cells selected for further analysis: ', ...
            inset,' of ',params.metric,' ',params.dimension];
        disp(msg)
    end
end