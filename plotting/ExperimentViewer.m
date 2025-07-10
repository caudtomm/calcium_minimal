classdef ExperimentViewer
    properties
        name
        subjectTab
        locations
        traces

        subjects_to_use logical
        filtered_traces

        dataFilter DataFilter
        plotConfig PlotConfig
    end

    methods
        function obj = ExperimentViewer(experiment)
            arguments
                experiment Experiment
            end
            obj.name = experiment.name;
            obj.subjectTab = experiment.subjectTab;
            obj.locations = experiment.locations;
            obj.traces = experiment.traces(:);
            obj.dataFilter = DataFilter();
            obj.plotConfig = PlotConfig();
        end

        %% getters

        function idx = get.subjects_to_use(obj)
            subjectIDs = obj.dataFilter.getSubjectIDs(obj.subjectTab);
            idx = ismember(obj.subjectTab.name,subjectIDs);
        end

        function traces = get.filtered_traces(obj)
            % filter subjects
            traces = obj.traces(obj.subjects_to_use);
        end

        %% setters

        function obj = setTheme(obj, themeName)
            obj.plotConfig.theme = themeName;
        end
        
        %% intermediate-level plotting function headers 
        % (normally call external low level functions that don't rely on custom objects)
        % still output a single plot onto provided axes

        function out = plotDistancesHead(obj, varargin)
            % plotDistances - Plot similarity/distance metrics for peri-stimulus traces across subjects.
            %
            % Usage:
            %   [hf, out] = obj.plotDistances('ps_lim', [start end], 'method', methodName, 'trial_sorting', sortingType)
            %
            % Inputs (as name-value pairs):
            %   'ps_lim'        - 2-element vector specifying peri-stimulus window in seconds [default: [1 20]]
            %   'plotType'      - String selecting which plot to produce
            %                     opitons: {'full','repetitions'}
            %   'method'        - String specifying similarity/distance metric (e.g., 'correlation') [default: 'correlation']
            %                     options: see pdist
            %   'trial_sorting' - String specifying trial sorting method (e.g., 'chronological') [default: 'chronological']
            %                     options: {'chronological','random','stim_id'}
            %   'stim_allowed'  - Cell array of strings containing labels of all stimuli to consider. (eg. {'Arg', 'Leu'})
            %                     or stimulus group name [default: 'all trials']
            %                     options: {'all trials','all stimuli','all CS+','all CS-','all familiar','all novel'}
            %
            % Outputs:
            %   out  - Output structure
            arguments
                obj
            end
            arguments (Repeating)
                varargin
            end

            % Set default values
            ps_lim = [1 20];
            plotType = 'full';
            method = 'correlation';
            trial_sorting = 'chronological';
            stim_allowed = 'all trials';

            % Parse name-value pairs
            if ~isempty(varargin)
                for k = 1:2:length(varargin)
                    switch lower(varargin{k})
                        case 'ps_lim'
                            ps_lim = varargin{k+1};
                        case 'plottype'
                            plotType = varargin{k+1};
                        case 'method'
                            method = varargin{k+1};
                        case 'trial_sorting'
                            trial_sorting = varargin{k+1};
                        case 'stims_allowed'
                            stim_allowed = varargin{k+1};
                    end
                end
            end

            % initialize output
            out = [];

            % Isolating relevant data
            nsubjects = numel(obj.filtered_traces);
            events = cell(nsubjects,1);
            all_labs = cell(nsubjects,1);
            for i = 1:nsubjects
                thistrace = obj.filtered_traces{i};

                % Trial sorting
                [~,trial_idx] = TraceViewer(thistrace).sortTrials(trial_sorting);

                % get peri-stimulus data [t,N,trials]
                M = thistrace.(obj.dataFilter.traceType)(:,:,trial_idx);
                stim_on_frame = thistrace.stim_series.frame_onset(1);
                fs = thistrace.framerate;
                events{i} = TraceViewer.getPeriEventData(M,stim_on_frame,ps_lim,fs);

                % Retrieve stimulus identity labels
                all_labs{i} = thistrace.stim_series.stimulus(trial_idx);

                % Stimulus filtering
                thisgroup = thistrace.subject_group;
                desired_stimuli = getStimuliByGroup(thisgroup,stim_allowed);
                idx = ismember(all_labs{i}, desired_stimuli);
                if ~isempty(desired_stimuli) && ~all(idx) && sum(idx)>0
                    % If some trials are not in the desired stimuli, filter them out
                    events{i} = events{i}(:,:,idx);
                    all_labs{i} = all_labs{i}(idx);
                elseif isempty(desired_stimuli) || sum(idx)==0
                    % If no stimuli are accepted, return without trying to
                    % plot ... nothing!
                    return
                end
            end

            % by default, labels are applied based on subject 1
            labs = all_labs{1};

            % Check that all label arrays are the same. This can happen
            % when filtering stimuli according to other properties than
            % identity (ex. trained/novel) from multiple experimental
            % groups.
            all_equal = all(cellfun(@(x) isequal(x, labs), all_labs));
            if ~all_equal
                warning('Mock labels! Actual stimulus identities differ across subjects!')
                % Replace labs with mock labels following the identity structure in labs{1}
                [~, ~, ic] = unique(labs, 'stable');
                mockLabels = arrayfun(@(x) char('A' + x - 1), ic, 'UniformOutput', false);
                labs = mockLabels;
            end

            % call low-level plotter
            out.distMat3d = plotDistances(events,plotType,method,labs,obj.plotConfig);
            title([num2str(ps_lim(1)),'-',num2str(ps_lim(2)), ' s'], ...
                'Color',obj.plotConfig.textcol)

            % return
            out.all_labs = all_labs;

        end


        function out = plotDiscriminationHead(obj, varargin)
            % plotDiscrimination - Plot classification performance / discriminability of a across subjects.
            %
            % Usage:
            %   [hf, out] = obj.plotDistances('ps_lim', [start end], 'method', methodName, 'trial_sorting', sortingType)
            %
            % Inputs (as name-value pairs):
            %   'ps_lim'        - 2-element vector specifying peri-stimulus window in seconds [default: [1 20]]
            %   'plotType'      - String selecting which plot to produce
            %                     opitons: {'full','repetitions'}
            %   'method'        - String specifying similarity/distance metric (e.g., 'correlation') [default: 'correlation']
            %                     options: see pdist
            %   'stim_allowed'  - Cell array of strings containing labels of all stimuli to consider. (eg. {'Arg', 'Leu'})
            %                     or stimulus group name [default: 'all trials']
            %                     options: {'all trials','all stimuli','all CS+','all CS-','all familiar','all novel'}
            %
            % Outputs:
            %   out  - Output structure
            arguments
                obj
            end
            arguments (Repeating)
                varargin
            end

            % Set default values
            ps_lim = [1 20];
            plotType = 'performance_lines';
            method = 'correlation';
            stim_allowed = 'all stimuli';
            reps_touse = [];
            focus_stims = 'all trials'; % by default, no further filtering
            do_zscore = false;

            % Parse name-value pairs
            if ~isempty(varargin)
                for k = 1:2:length(varargin)
                    switch lower(varargin{k})
                        case 'ps_lim'
                            ps_lim = varargin{k+1};
                        case 'plottype'
                            plotType = varargin{k+1};
                        case 'method'
                            method = varargin{k+1};
                        case 'repetitions'
                            reps_touse = varargin{k+1};
                        case 'stims_allowed'
                            stim_allowed = varargin{k+1};
                        case 'focus_stims'
                            focus_stims = varargin{k+1};
                        case 'zscore'
                            do_zscore = varargin{k+1};
                    end
                end
            end

            % initialize output
            out = [];

            % Isolating relevant data
            nsubjects = numel(obj.filtered_traces);
            events = cell(nsubjects,1);
            all_labs = cell(nsubjects,1);
            for i = 1:nsubjects
                thistrace = obj.filtered_traces{i};

                % get peri-stimulus data [t,N,trials]
                M = thistrace.(obj.dataFilter.traceType);
                stim_on_frame = thistrace.stim_series.frame_onset(1);
                fs = thistrace.framerate;
                events{i} = TraceViewer.getPeriEventData(M,stim_on_frame,ps_lim,fs);

                % Retrieve stimulus identity labels
                all_labs{i} = thistrace.stim_series.stimulus;

                % Stimulus filtering
                thisgroup = thistrace.subject_group;
                desired_stimuli = getStimuliByGroup(thisgroup,stim_allowed);
                idx = ismember(all_labs{i}, desired_stimuli);
                if ~isempty(desired_stimuli) && ~all(idx)
                    % If some trials are not in the desired stimuli, filter them out
                    events{i} = events{i}(:,:,idx);
                    all_labs{i} = all_labs{i}(idx);
                elseif isempty(desired_stimuli) || sum(idx)==0
                    % If no stimuli are accepted, return without trying to
                    % plot ... nothing!
                    return
                end

                % Stimulus repetition filter
                if isempty(reps_touse); continue; end % empty argument 'repetitions' leads to all repetitions being used
                thisstims = unique(all_labs{i});
                nstims = numel(thisstims);
                idx_keep = false(1, numel(all_labs{i}));
                for i_stim = 1:nstims
                    idx_stim = find(ismember(all_labs{i}, thisstims{i_stim}));
                    this_nreps = numel(idx_stim);
                    % Select only allowed repetition indices
                    reps_available = 1:this_nreps;
                    reps_valid = reps_available(ismember(reps_available, reps_touse));
                    if isempty(reps_valid)
                        continue
                    end
                    idx_keep(idx_stim(reps_valid)) = true;
                end
                if ~all(idx_keep)
                    % Filter events and labels to keep only desired repetitions
                    events{i} = events{i}(:,:,idx_keep);
                    all_labs{i} = all_labs{i}(idx_keep);
                elseif sum(idx_keep)==0
                    % If no trials are accepted, return without trying to
                    % plot ... nothing!
                    return
                end
            end

            % call low-level processor (perform discrimination analysis)
            all_out = cell(nsubjects,1);
            for i = 1:nsubjects
                thisevents = events{i};
                thislabs = all_labs{i};

                % if there are no allowed stimuli here, skip subject
                if isempty(thislabs); continue; end

                % call post-processing function
                all_out{i} = doDiscrimination(thisevents, thislabs, ...
                                                       'method', method);
            end

            % focus on specific trials for plotting (without changing any of the values!)
            focus_trials = cell(nsubjects,1);
            for i = 1:nsubjects
                thistrace = obj.filtered_traces{i};
                thisgroup = thistrace.subject_group;
                desired_stimuli = getStimuliByGroup(thisgroup,focus_stims);
                idx = ismember(all_labs{i}, desired_stimuli);
                focus_trials{i} = find(idx);
            end

            % get rid of empty data
            idx = cellfun(@isempty,all_out) | cellfun(@isempty,focus_trials);
            all_out(idx) = [];
            focus_trials(idx) = [];
            if isempty(all_out); return; end

            % call low-level plotter
            out = plotDiscrimination(all_out,...
                plotType,obj.plotConfig, ...
                'method',method, ...
                'FocusTrials',focus_trials, ...
                'actualreps', reps_touse, ...
                'zscore',do_zscore);

            % return

        end



        %% complex and idiosyncratic high-level plotters

        function [hf, out] = plotRepetitionDistances(obj, ps_lim, method)  
            arguments
                obj
                ps_lim = [1 20]
                method = 'correlation'
            end
            
            % Knobs
            groups = {'naïve', ...
                      'trained', ...
                      'trained1', ...
                      'trained2', ...
                      'trained1 (T-R-S-H-A-ACSF/L)', ...
                      'uncoupled'};
            odor_sets = {'all stimuli', ...
                         {'Arg','Ala','His','Trp','Ser'}, ...
                         {'Leu'}, ...
                         'all CS+', ...
                         'all CS-', ...
                         'all familiar', ...
                         'all novel', ...
                        };

            % useful metrics
            ngroups = numel(groups);
            nodor_sets = numel(odor_sets);
            nplots = ngroups*nodor_sets;

            % Initialize output
            hf = gobjects(2,1);
            out = cell(nplots,1);


            %% Figure 1: comparison of distance matrices

            % define figure size
            ncols = nodor_sets;
            nrows = ngroups;

            hf(1) = figure; % [groups, odor_sets]
            n = 1;
            labs = cell(nplots,1);
            for i_g = 1:ngroups
                % filter data by group
                thisgroup = groups{i_g};
                obj.dataFilter.subjectGroup = thisgroup;

                for i_o = 1:nodor_sets
                    thisodorset = odor_sets{i_o};

                    thisodorset_str = thisodorset;
                    if iscell(thisodorset_str); thisodorset_str = strjoin(thisodorset, ', '); end
                    msg = ['Plotting group ''',thisgroup,''' for odors: ',thisodorset_str];
                    disp(msg)

                    % build axes
                    subplot(nrows,ncols,n);

                    % call intermediate-level plotter
                    out{n} = obj.plotDistancesHead('ps_lim',ps_lim, ...
                                    'plotType', 'repetitions', ...
                                    'method',method, ...
                                    'stims_allowed',thisodorset);

                    % override title and ylabel
                    title(thisodorset_str)
                    ylabel(thisgroup)

                    % if there is no data, delete the subplot
                    if isempty(out{n}); axis off; end

                    % export label
                    labs{n} = [thisgroup,' - ',thisodorset_str];

                    % advance axis counter
                    n = n+1;
                end
            end
            set(gcf,'Position',[1 1 2000 1000])

            %% Figure 2: comparison of distance distributions

            % compute data
            data = [];
            for i = 1:nplots
                % skip empty output
                if isempty(out{i}); continue; end

                thismat = 1 - out{i}.distMat3d; % convert distance to similarity
                
                % knobs # TODO : tunable param
                repetitions = 2:5;

                % crop and take out the diagonal
                thismat = thismat(repetitions,repetitions,:);
                identity = logical(eye(length(repetitions)));
                thismat(identity) = nan;

                % store
                plot_idx = i * ones(numel(thismat),1);
                data = [data; plot_idx, thismat(:)];
            end
            
            % plot
            cfg = obj.plotConfig;
            hf(2) = figure;
            boxplot(data(:,2),labs(data(:,1)),'Orientation','horizontal','PlotStyle','compact','Colors',cfg.textcol)
            box off
            set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol);
            set(gcf, 'color', cfg.bgcol); 
            xlabel(method,'Color',cfg.textcol)
            set(gcf,'Position',[1 1 2000 1000])
            hold off
        end

        function [hf, out] = plotDiscriminationPerformanceMats(obj, ps_lim, method,focus_stims,repetitions,do_zscore)  
            arguments
                obj
                ps_lim = [1 20]
                method = 'correlation'
                focus_stims = 'all trials' % by default no further filtering
                repetitions = [];
                do_zscore = false
            end
            
            % Knobs
            groups = {'naïve', ...
                      'trained', ...
                      'trained1', ...
                      'trained2', ...
                      'trained1 (T-R-S-H-A-ACSF/L)', ...
                      'uncoupled'};
            odor_sets = {'all stimuli', ...
                         {'Arg','Ala','His','Trp','Ser'}, ...
                         'all familiar', ...
                         'all novel', ...
                        };

            % useful metrics
            ngroups = numel(groups);
            nodor_sets = numel(odor_sets);
            nplots = ngroups*nodor_sets;

            % Initialize output
            hf = gobjects(1,1);
            out = cell(nplots,1);


            %% Figure 1: comparison of discrimination performance curves

            % define figure size
            ncols = nodor_sets;
            nrows = ngroups;

            hf(1) = figure; % [groups, odor_sets]
            n = 1;
            labs = cell(nplots,1);
            for i_g = 1:ngroups
                % filter data by group
                thisgroup = groups{i_g};
                obj.dataFilter.subjectGroup = thisgroup;

                for i_o = 1:nodor_sets
                    thisodorset = odor_sets{i_o};

                    thisodorset_str = thisodorset;
                    if iscell(thisodorset_str); thisodorset_str = strjoin(thisodorset, ', '); end
                    msg = ['Plotting group ''',thisgroup,''' for odors: ',thisodorset_str];
                    disp(msg)

                    % build axes
                    subplot(nrows,ncols,n);

                    % call intermediate-level plotter
                    out{n} = obj.plotDiscriminationHead('ps_lim',ps_lim, ...
                                    'plotType', 'performance_mat', ...
                                    'method',method, ...
                                    'focus_stims', focus_stims, ...
                                    'zscore',do_zscore, ...
                                    'repetitions',repetitions, ...
                                    'stims_allowed',thisodorset);

                    % override title and ylabel
                    title(thisodorset_str,'color',obj.plotConfig.textcol)
                    ylabel(thisgroup)

                    % if there is no data, delete the subplot
                    if isempty(out{n}); axis off; end

                    % export label
                    labs{n} = [thisgroup,' - ',thisodorset_str];

                    % advance axis counter
                    n = n+1;
                end
            end
            set(gcf,'Position',[1 1 2000 1000])


        end
        
        % compare dFoverF and pSpike
        function [hf, cell_metrics] = compareTraceTypes(obj)
            % knobs
            nbins = 30; % for deviation from linearity (applied to pSpike)
                        % and avg transfer function curve (applied to dF)
            ncellscatterplots = 3;
            sort_bygroup = true;
            
            % init vars
            nfish = numel(obj.traces);
            hf = gobjects(9,1);
            
            % figure init
            hf(1) = figure;
            ncols = floor(sqrt(nfish))+1;
            nrows = ncols-1;
            
            cell_metrics = cell(nfish, 1);
            dF_all = [];
            pSpike_all = [];
            for i = 1:nfish
                disp(['fish #',num2str(i)])
                
                traces = obj.traces{i};
                dF = traces.format(traces.dFoverF_good);
                pSpike = traces.format(traces.pSpike);
            
                [T,N] = size(dF);
                disp([num2str(T),' timepoints, ',num2str(N),' cells.'])
                
                % pSpike vs dF for each cell 
                cells_toplot = randi(N,ncellscatterplots,1);
                polycoef = zeros(N,2);      % linear fit coefficients
                deviations = zeros(N,1);    % 'deviation from linearity' (see Rupprecht et al. 2025)
                ndatapoints = zeros(N,1);
                normr = zeros(N,1); 
                curves = zeros(nbins,N);
                for i_cell = 1:N
                    thisdF = dF(:,i_cell);
                    thispSpike = pSpike(:,i_cell);
            
                    % clean up nan values
                    to_rem = isnan(thisdF) | isnan(thispSpike);
                    thisdF = thisdF(~to_rem);
                    thispSpike = thispSpike(~to_rem);
                    ndatapoints(i_cell) = length(thisdF);
                    dF_all = [dF_all; thisdF];
                    pSpike_all = [pSpike_all; thispSpike];
            
                    % linear fit
                    [p,S] = polyfit(thisdF, thispSpike, 1);
                    polycoef(i_cell,:) = p;
            
                    % (debug) plot single cell scatter
                    % if ismember(i_cell,cells_toplot)
                    %     figure; scatter(thisdF,thispSpike,4,'k','filled')
                    %     axis square tight;
                    %     x = xlim; y = polyval(p,x);
                    %     hold on; line(x,y,'Color','r','LineWidth',1.5)
                    %     xlabel('dFoverF'); ylabel('inferred SR')
                    %     title(['cell #',num2str(i_cell)])
                    % end
                    
                    % deviation from linearity
                    thispSpike_linear = arrayfun(@(x) p(1)*x+p(2),thisdF);
                    bin_edges = linspace(min(thispSpike),max(thispSpike),nbins+1);
                    deviations_bins = zeros(nbins,1);
                    for i_bin = 1:nbins
                        % get average values in the bin
                        idx = thispSpike>=bin_edges(i_bin) & thispSpike<=bin_edges(i_bin+1);
                        binpSpike = mean(thispSpike(idx));
                        binpSpike_linear = mean(thispSpike_linear(idx));
                        
                        % get deviation
                        deviations_bins(i_bin) = (binpSpike - binpSpike_linear)^2 / (binpSpike_linear^2);
                    end
                    deviations(i_cell) = sqrt(mean(deviations_bins,'omitmissing'));
            
                    % norm of fit residuals
                    normr(i_cell) = S.normr;
            
                    % avg transfer function curve
                    bin_edges = linspace(min(dF,[],'all'),max(dF,[],'all'),nbins+1); % bins over all cells
                    for i_bin = 1:nbins
                        % get average value in the bin
                        idx = thisdF>=bin_edges(i_bin) & thisdF<=bin_edges(i_bin+1);
                        curves(i_bin,i_cell) = mean(thispSpike(idx));
                    end
                end
            
                % store single-cell metrics to table
                slopes = polycoef(:,1);
                intercepts = polycoef(:,2);
                noise_level_PR = traces.dFnoise(traces.goodNeuron_IDs);
                pxVariance_overtime = traces.format(traces.Fnoise(:,traces.goodNeuron_IDs,:));
                pxVariance = mean(pxVariance_overtime,'omitmissing');
                pxVariance = pxVariance(:);
                cell_metrics{i} = table( ...
                    ndatapoints, ...
                    deviations, ...
                    slopes, ...
                    intercepts, ...
                    normr, ...
                    noise_level_PR, ...
                    pxVariance);
            
                % transfer function plot
                figure(hf(1));
                subplot(nrows,ncols,i)
                x = repmat(bin_edges(2:end),N,1)'; % [nbins,N]
                plot(x,curves,'k')
                axis square tight; box off
                xlabel('dFoverF'); ylabel('inferred SR');
                title(['fish #',num2str(i)])
            end
            
            % sort order of fish by group
            [sorted_groups,idxbygroup] = sort(obj.subjectTab.group);
            cell_metrics_sorted = cell_metrics(idxbygroup);
            grouplabels = obj.subjectTab.group;
            
            % cell metrics boxplots
            labels = cell_metrics{1}.Properties.VariableNames;
            for i = 1:numel(labels)
                thislabel = labels{i};
                values = [];
                for i_fish = 1:nfish
                    thisvalues = cell_metrics_sorted{i_fish}.(thislabel);
                    N = numel(thisvalues);
                    values = [values; thisvalues repelem(i_fish,N,1)];
                end
                hf(1+i) = figure;
                boxplot(values(:,1),values(:,2),'PlotStyle','compact','Symbol','');
                xticks(1:nfish); xticklabels(grouplabels);xtickangle(90)
                title(thislabel)
                box off
            end
            
            % heatmap of inferred SR vs dFoverF over all data
            hf(9) = figure;
            subplot(131); h(1)=histogram(dF_all,50);
            xlabel('dFoverF'); axis square; box off
            subplot(132); h(2)=histogram(pSpike_all,50);
            xlabel('inferred SR'); axis square; box off
            x1 = h(1).BinEdges(2:end); x2 = h(2).BinEdges(2:end);
            vals1 = h(1).Values; vals2 = h(2).Values;
            hmap = vals2' * vals1;
            subplot(133); imagesc(x1,x2,log(hmap));
            axis square
            b = colorbar; b.Label.String = 'Log Density';
            xlabel('dFoverF'); ylabel('inferred SR');


        end

    end
end
