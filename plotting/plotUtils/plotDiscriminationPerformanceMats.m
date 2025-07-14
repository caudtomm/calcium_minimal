
        function [hf, out, groups, odor_sets] = plotDiscriminationPerformanceMats(v, ps_lim, method,focus_stims,repetitions,do_zscore)  
            arguments
                v ExperimentViewer
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
                v.dataFilter.subjectGroup = thisgroup;

                for i_o = 1:nodor_sets
                    thisodorset = odor_sets{i_o};

                    thisodorset_str = thisodorset;
                    if iscell(thisodorset_str); thisodorset_str = strjoin(thisodorset, ', '); end
                    msg = ['Plotting group ''',thisgroup,''' for odors: ',thisodorset_str];
                    disp(msg)

                    % build axes
                    subplot(nrows,ncols,n);

                    % call intermediate-level plotter
                    out{n} = v.plotDiscriminationHead('ps_lim',ps_lim, ...
                                    'plotType', 'performance_mat', ...
                                    'method',method, ...
                                    'focus_stims', focus_stims, ...
                                    'zscore',do_zscore, ...
                                    'repetitions',repetitions, ...
                                    'stims_allowed',thisodorset);

                    % override title and ylabel
                    title(thisodorset_str,'color',v.plotConfig.textcol)
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
        