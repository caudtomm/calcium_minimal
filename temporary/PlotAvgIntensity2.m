%% set run parameters
run_grid = false;

%% single run parameters
params = Params(); % update attributes below
params.groups = {'naïve';'trained1';'trained2';'uncoupled'}; ... % options : {'previousnaive';'naïve';'trained1';'trained2';'trained1alt';'uncoupled'};
params.stimuli = {'Arg','Ala','His','Trp','Ser','Leu','ACSF','spont.'}; ... % options : {'all CS+','all CS-','all familiar','all novel'} or {'Arg','Ala','His','Trp','Ser','Leu','Food','ACSF','spont.'};
params.time_range = [1 20]; ... % time range : [from stim_on, from stim_on] (sec)
params.only_top_variant_units = false; ...  % only_top_variant_units
params.metric = 'avg_activity'; ... % metric : options {'max_activity','avg_activity','variance'}
params.dimension = 'over_time'; ...% options: {'over_cells', 'over_time'}
params.n_equals = 'all'; ... % n = ?
params.order_by = 'trial_num'; ... % order by : options {'trial_num','relative_trial_num','stim_type'};
params.todo_reltrialnum = [1:5]; % todo_reltrialnum
params.isolate_quantile = 'all';
params.quantile = .9;

%% 'all_options' structure
all_options.groups = {...
    {'previousnaive';'naïve';'trained1';'trained2';'trained1alt'}};
all_options.stimuli = {...
    'all CS+', ...
    'all CS-', ...
    'all familiar', ...
    'all novel', ...
    {'Arg','Ala','His','Trp','Ser','Leu','Food'}, ...
    {'Arg','Ala','His','Trp','Ser','Leu','Food','ACSF','spont.'}, ...
    'Leu'};
all_options.time_range = {...
    [0,20], ...
    [0,3],...
    [3,9], ...
    [12,20], ...
    [21,25]};
all_options.only_top_variant_units = {true, false};
all_options.metric = {'max_activity','avg_activity','variance'};
all_options.dimension = {'over_cells', 'over_time'};
all_options.n_equals = {'fish'};
all_options.order_by = {'relative_trial_num'};
all_options.y_range = {[0 5]};
all_options.todo_reltrialnum = {[1:5],[1:12]};


%% run

if run_grid
    runGridSearch(experiment, all_options);
else
    [datatab,h] = plotActivityMetric(experiment, params);
end


%% Main functions

function [datatab,ax] = plotActivityMetric(experiment,params)
%% get plottable data table [organizer,metric]
[output_struct, datatab] = extractActivityMetric(experiment,params);
n_vals = size(datatab,1);

%% plot the progression of the metric across relative trials

figlab = [strcat(params.stimuli),', ',params.metric,' ',params.dimension];

% take average plot line
orgs = sort(unique(datatab.organizer_value));
n_orgvals = numel(orgs);
ymean = nan(numel(orgs),1);
ymedian = nan(numel(orgs),1);
for i = 1:n_orgvals
    idx = ismember(datatab.organizer_value,orgs(i));
    ymean(i) = nanmean(datatab.metric_value(idx));
    ymedian(i) = nanmedian(datatab.metric_value(idx));
end
ymean_orgtrue = nan(max(orgs),1);
ymedian_orgtrue = nan(max(orgs),1);
for i=1:n_orgvals
    ymean_orgtrue(orgs(i)) = ymean(i);
    ymedian_orgtrue(orgs(i)) = ymedian(i);
end

% exponential fit
% !!!! starting nan because first trial is baseline!
% n_exclude = 7
% fit_to = [ nan(1,min([max(orgs),n_exclude])) , ...
%     ymean_orgtrue(min([max(orgs),n_exclude+1]):min([max(orgs),24]))' ];
% samplerange = 1:max(orgs)-1;
% [efit] = [nan,fitExponentialDecay(fit_to,samplerange)];
% % [efit] = [nan,fitExponentialDecay([nan(1,n_exclude),ymean_orgtrue(n_exclude+1:end)'],1:max(orgs)-1)];
% 
% % add data column: metric - exp. fit
% tmp = nan(n_vals,1);
% for i=1:n_vals
%     tmp(i) = datatab.metric_value(i) - efit(1,datatab.organizer_value(i));
% end
% datatab.metric_minus_efit = tmp;

ax= [];
% linear plot
% figure
% scatter(datatab,"organizer_value","metric_value")
figure;
boxplot(datatab.metric_value,datatab.organizer_value,'Color','w','symbol' ,' ')
hold on; scatter(1:numel(ymean),ymean,80,'r','filled')
% plot(efit(~isnan(ymedian_orgtrue)),'r-','LineWidth',2);
hold on
title(figlab,'Color','k')
try; ylim([0 1.2*nanmax(ymean)]); catch; end
set(gca, 'color', 'none', 'XColor','w', 'YColor','w', 'ZColor','w','FontSize',14);
set(gcf, 'color', 'none'); 
set(gcf, 'Position', [50 50 350 400]);
ax = gca;
box off

%loglog plot
figure;
loglog(1:numel(ymean_orgtrue),ymean_orgtrue,'Color','w');hold on
% loglog(efit,'r-','LineWidth',2); hold on
title(figlab,'Color','w')
try; ylim([0 1.2*nanmax(ymean_orgtrue)]); catch; end
set(gca, 'color', 'none', 'XColor','w', 'YColor','w', 'ZColor','w','FontSize',14);
set(gcf, 'color', 'none'); 
set(gcf, 'Position', [50 50 350 400]);

figure;
% plot(efit,'r-','LineWidth',2); hold on
err = nan(1,numel(ymean_orgtrue));
for i=1:max(orgs)
    idx = datatab.organizer_value==i;
    if isempty(idx); continue; end
    tmp = datatab(idx,:);
    err(i)=std(tmp.metric_value);
end
errorbar(ymean_orgtrue,err,'Color','k');
title(figlab,'Color','k')
try; ylim([0 1.2*nanmax(ymean)]); catch; end
xlim([.5, .5+numel(ymean)])
set(gcf, 'Position', [50 50 350 400]);

figure;
plot(zeros(1,max(orgs)),'b--','LineWidth',3); hold on
% scatter(datatab,"organizer_value","metric_minus_efit",'filled','SizeData',5,'MarkerFaceColor','cyan')
% plot(ymean_orgtrue'-efit,'r-','LineWidth',5)
% ylim([-3 25])
xlabel(params.order_by)
ylabel([params.metric,' - fit'])



end

function runGridSearch(experiment, all_options)
    % Extract field names and initialize a cell array to hold the option combinations
    options = combine_options(all_options);
    
    % Assuming options is a table where each row represents a parameter combination
    num_param_sets = height(options);
    
    % Grid size, page size, and filename as before
    gridSize = [3, 4];
    pageSize = [8.27, 11.69]; % A4 paper size in inches
    filename = 'compiled_figures.pdf';
    
    % Figures per page
    figures_per_page = prod(gridSize);
    h = [];
    
    % Loop over all parameter combinations
    for p = 1:num_param_sets
        
        % Extract parameters for this combination
        params = storeOptions2Params(options(p, :));
        
        [~,ax] = plotActivityMetric(experiment, params);
        h{p} = ax(1);
        
        % Generate a title string from the parameter combination
        param_str = generate_param_str(params);
        title(['Parameters: ' param_str(1:end-2)]);
    
    end
    
    num_figures = numel(h);
        
    % Compilation of figures into PDF, with adjusted space for readability
    for i = 1:ceil(num_figures/figures_per_page)
        
        % Create a new figure for the compiled plots, invisible, with A4 size
        h_compiled = figure('Visible', 'off', 'PaperPosition', [0 0 pageSize], ...
                            'PaperSize', pageSize, 'PaperUnits', 'inches');
        % Ensure h_compiled is the current figure
        figure(h_compiled);
    
        for j = 1:figures_per_page
            idx = (i-1)*figures_per_page + j;
            if idx <= num_figures
                % Ensure h_compiled is the current figure
                figure(h_compiled);
                
                % Create subplot, ensure ax is an Axes object.
                ax = subplot(gridSize(1), gridSize(2), j, 'Parent', h_compiled);
                
                % Double-check object to copy is an axes, and if so, copy children.
                if isaxes(h{idx})  % Using isaxes to ensure object is an axes
                    copyobj(allchild(h{idx}), ax);
                else
                    warning('Object #%d is not an axes and was not copied.', idx);
                end
                
                % Remove title (since params are summarized in PDF name)
                title(ax, num2str(idx));
                
                % Optionally: adjust ax.Position here for more space, if needed
            end
        end
        set(gcf, 'color', 'w');
    
        % Export to PDF, appending if file already exists
        if exist(filename, 'file')
            export_fig(filename, h_compiled, '-pdf', '-nocrop', '-append');
        else
            export_fig(filename, h_compiled, '-pdf', '-nocrop');
        end
        
        % Close the compiled figure
        close(h_compiled);
        
    end
       
    save('compiled_figures_options.mat', 'options');
end


%% Accessory functions

% Recursive function to generate all combinations of options
function options_table = combine_options(all_options)
% combine_options: Generates a table of all unique combinations of options.
%
% INPUT:
%    all_options: a structure, each field of which is a cell array of options.
%
% OUTPUT:
%    options_table: a table with one column for each field of 'all_options' 
%                   and one row for each unique combination of field options.

    % Extract field names and number
    fields = fieldnames(all_options);
    num_fields = numel(fields);
    all_combinations = cell(0, num_fields);

    % Recursive function to generate all combinations of options
    function combine_options_recursive(prefix, field_idx)
        if field_idx > num_fields
            all_combinations(end+1, :) = prefix;
            return;
        end

        options = all_options.(fields{field_idx});
        for i = 1:numel(options)
            combine_options_recursive([prefix, options(i)], field_idx + 1);
        end
    end

    % Initialize the recursive function
    combine_options_recursive({}, 1);

    % Convert the cell array to a table
    options_table = cell2table(all_combinations, 'VariableNames', fields);
end

% Function to construct a string from parameter names and values.
function param_str = generate_param_str(params)
% generate_param_str: Constructs a string from parameter names and values.
%
% INPUT:
%    params: a one-row table containing parameters and their values.
%
% OUTPUT:
%    param_str: a string, containing parameter names and their corresponding values.

    % Initialize an empty string
    param_str = '';

    % Loop over each variable in the table, adding each name-value pair to the string
    for k = 1:width(params)
        var_name = params.Properties.VariableNames{k};
        var_val = params{1, k};

        % Convert the variable value to a string based on its type
        if iscell(var_val); var_val=var_val{1}; end
        if isnumeric(var_val)
            val_str = mat2str(var_val);
        elseif iscell(var_val)
            if all(cellfun(@isnumeric, var_val))
                val_str = ['{' strjoin(cellfun(@mat2str, var_val, 'UniformOutput', false), ', ') '}'];
            elseif all(cellfun(@(x) isstring(x) || ischar(x), var_val))
                val_str = ['{' strjoin(cellfun(@char, var_val, 'UniformOutput', false), ', ') '}'];
            elseif all(cellfun(@islogical, var_val))
                val_str = ['{' strjoin(cellfun(@logical2str, var_val, 'UniformOutput', false), ', ') '}'];
            else
                error('Unsupported cell array content.');
            end
        elseif isstring(var_val)
            val_str = char(var_val);
        elseif ischar(var_val)
            val_str = var_val;
        elseif islogical(var_val)
            val_str = logical2str(var_val);
        else
            error('Unsupported data type.');
        end

        % Add the name-value pair to the string
        param_str = [param_str, sprintf('%s: %s, ', var_name, val_str)];
    end

    % Remove the trailing comma and space
    param_str = param_str(1:end-2);
end

function str = logical2str(log_val)
    % Convert a logical value to a string ('true' or 'false')
    if log_val
        str = 'true';
    else
        str = 'false';
    end
end

function isAx = isaxes(h)
    % isaxes: Checks whether the input handle(s) h are axes.
    %
    % INPUT:
    %    h: An array of graphic handle(s).
    %
    % OUTPUT:
    %    isAx: Logical array, true where handle corresponds to an axes.

    % Initialize logical array
    isAx = false(size(h));
    
    % Check if each handle is an axes
    for k = 1:numel(h)
        isAx(k) = isgraphics(h(k), 'axes');
    end
end

function params = storeOptions2Params(options)
    %assumes single-line table input
    params.groups = options.groups{1};
    params.stimuli = options.stimuli{1};
    params.time_range = options.time_range;
    params.only_top_variant_units = options.only_top_variant_units;
    params.metric = options.metric{1};
    params.dimension = options.dimension{1};
    params.n_equals = options.n_equals{1};
    params.order_by = options.order_by{1};
    params.y_range =  options.y_range;
    params.todo_reltrialnum = options.todo_reltrialnum{1};
end


