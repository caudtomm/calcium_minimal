function [hf, out, groups, odor_sets] = plotModeDecompositionFigure(v, ps_lim, method, plotType, varargin)  
    arguments
        v ExperimentViewer
        ps_lim = [1 20]
        method = 'nmf'
        plotType = 'curves' % 'curves' or 'trialmap'
    end
    
    % Knobs
    groups = {'naïve', ...
                'trained', ...
                'uncoupled'};
    odor_sets = {'all stimuli', ...
                'all familiar', ...
                'all novel', ...
                };

    % useful metrics
    ngroups = numel(groups);
    nodor_sets = numel(odor_sets);
    nplots = ngroups*nodor_sets;

    % Initialize output
    hf = gobjects(10,1);
    i_hf = 0;
    out = cell(nplots,1);

    %% Figure 1: low dimensional embedding of single unit tuning across repetitions,
    % based on common axes ( one figure for each subplot from Figure 1)
    
    for i_g = 1:ngroups
        thisgroup = groups{i_g};
        v.dataFilter.subjectGroup = thisgroup;
        
        for i_os = 1:nodor_sets
            thisodorset = odor_sets{i_os};
            
            i_hf = i_hf+1;
            hf(i_hf) = figure;
            cfg = v.plotConfig; % # TODO: use cfg to set up the figure

            out = v.getModeDecomposition('ps_lim', [1 20], ...
                                    method, 'nmf', ...
                                    trial_sorting, 'chronological', ...
                                    reps_touse, [], ...
                                    stim_allowed, 'all stimuli');
            if isempty(out); continue; end

            % Plot the results
            switch lower(plotType)
                case 'curves'
                    plotCurves();
                case 'trialmap'
                    plotTrialMap();
                otherwise
                    error('Unknown plot type: %s', plotType);
            end

            title([thisgroup,', ',thisodorset])
        end
    end

    function plotCurves()
        % Plot curves for the current group and odor set
        % This function should implement the logic to plot the curves
        % based on the output from doModeDecomposition.
        % Example:
        % plot(out.vals);
        xlabel('Time');
        ylabel('Activity');
        legend('show');
    end

    function plotTrialMap()
        % Plot trial map for the current group and odor set
        % This function should implement the logic to plot the trial map
        % based on the output from doModeDecomposition.
        % Example:
        % imagesc(out.vals);
        colorbar;
        xlabel('Trials');
        ylabel('Activity');
    end

end
