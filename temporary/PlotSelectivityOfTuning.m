
todo_group = {'naive'};

method = 'population';
% REMINDER: for 'kurtosis-stimuli', you need each fish to have the same
% number of odors to compare - EDIT 'stims2use' accordingly below

%%
todo_fish = [];
for i_fish = 1:numel(experiment.series)
    if ismember(experiment.series{i_fish}.group,todo_group)
        todo_fish = [todo_fish;i_fish];
    end
end


% ... and sparseness matrix [trials x fish]
population_sparsenessCELL = {};

%% Calculate and store sparseness for each trial and fish

for i_fish = 1:numel(todo_fish)
    thisgroup = experiment.series{todo_fish(i_fish)}.group;
    data = experiment.series{todo_fish(i_fish)}.data;
    fs = data.meta.framerate;

    switch thisgroup
        case 'previousnaive'
            stims2use = {'Trp','Ala','Ser','Food'};
        case {'trained1','naive','trained2','uncoupled'}
            stims2use = {'Arg','Ala','His','Trp'};
%             stims2use = {'Arg','Ala','His'};
%             stims2use = {'Trp','Ser','Leu'};
    end
    interval = [floor(data.stim_on_sec*fs) : 1+floor(data.stim_off_sec*fs)];

    % extract activity traces (only for odor periods)
    activityTraces = traceFormat(selectTimeFrame(data.traces,interval,data.L),length(interval));
    
    % eliminate any baseline or blank trials, and trials of unwanted
    % stimuli
    touse_trial = ismember(data.stim_type, stims2use);
    activityTraces(:,:,~touse_trial) = [];
    
    % calculate tuning selectivity
    sparseness = calculateSparseness(activityTraces,method,'StimTypes',data.stim_type(touse_trial));
    

    % store sparseness in output matrix
    population_sparsenessCELL{i_fish} = sparseness; 
end

% data for plotting
dataCELL = population_sparsenessCELL;
% labs = reference_data.stim_type(touse_trial);

%% plot histogram over fish

data = dataCELL;


nbins = 20;

range = [0,1];
if strcmp(method,'kurtosis-stimuli') | strcmp(method,'kurtosis-singletrials')
    range = [min(cellfun(@min,data))-1,max(cellfun(@max,data))+1];
end
edges = linspace(range(1),range(2),nbins+1);

% Create a matrix to store the item counts for each bin
itemCounts = zeros(numel(data), nbins);

% Calculate item counts for each dataset
for i_fish = 1:numel(data)
    counts = histcounts(data{i_fish},edges);
    itemCounts(i_fish, 1:numel(counts)) = counts;
end

% Normalize the item counts by the size of the vector
normalizedCounts = itemCounts ./ sum(itemCounts, 2);


figure;

% Plot the normalized counts
plot(edges(2:end),normalizedCounts','LineWidth',2);
hold on;

% Set the axis labels and title
xlabel(method);
ylabel('Norm. Histogram');

% Optionally, you can adjust the figure properties
set(gcf, 'Color', 'k');  % Set black background
set(gca, 'color', 'none', 'XColor','w', 'YColor','w', 'ZColor','w', 'FontSize',12);
set(gcf, 'Position', [50 50 300 250]);
axis tight
box off
hold off;


figure;

% Plot the normalized counts
semilogy(edges(2:end),normalizedCounts','LineWidth',2);
hold on;

% Set the axis labels and title
xlabel(method);
ylabel('Norm. Histogram');

% Optionally, you can adjust the figure properties
set(gcf, 'Color', 'k');  % Set black background
set(gca, 'color', 'none', 'XColor','w', 'YColor','w', 'ZColor','w');
set(gcf, 'Position', [50 50 400 350]);
axis tight
box off
hold off;

