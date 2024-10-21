classdef Params
    properties
        groups
        stimuli
        time_range
        only_top_variant_units
        metric
        dimension
        n_equals
        order_by
        todo_reltrialnum
        isolate_quantile
        quantile {mustBeInRange(quantile,0,1)}
    end
    methods (Static)
        % Constructor
        function params = Params()
            params.groups = {'previousnaive';'naive';'trained1';'trained2';'trained1alt'}; ... % options : {'previousnaive';'naive';'trained1';'trained2';'trained1alt';'trainedUC'};
            params.stimuli = {'Arg','Ala','His','Trp','Ser','Leu','Food'}; ... % options : {'all CS+','all CS-','all familiar','all novel'} or {'Arg','Ala','His','Trp','Ser','Leu','Food','ACSF','spont.'};
            params.time_range = [0 20]; ... % time range : [from stim_on, from stim_on] (sec)
            params.only_top_variant_units = false; ...  % only_top_variant_units
            params.metric = 'variance'; ... % metric : options {'max_activity','avg_activity','variance'}
            params.dimension = 'over_time'; ...% options: {'over_cells', 'over_time'}
            params.n_equals = 'fish'; ... % n = ?
            params.order_by = 'relative_trial_num'; ... % order by : options {'trial_num','relative_trial_num','stim_type'};
            params.todo_reltrialnum = [1:5]; % todo_reltrialnum
            params.isolate_quantile = 'all'; % options: {'all', 'bottom', 'top'}
            params.quantile = .9;
        end
    end
end