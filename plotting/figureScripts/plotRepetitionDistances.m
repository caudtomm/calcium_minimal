
function [hf, out, subject_groups, stim_groups] = plotRepetitionDistances(v, method, subject_groups, stim_groups, repetitions)
    arguments
        v ExperimentViewer
        method string = 'correlation'
        subject_groups cell = {'naïve', 'trained'}
        stim_groups cell = {{'Arg','Ala','His'}, {'Trp','Ser','Leu'}, 'all stimuli'}
        repetitions double = 1:5
    end

    ngroups    = numel(subject_groups);
    nodor_sets = numel(stim_groups);
    nplots     = ngroups * nodor_sets;

    hf  = gobjects(2,1);
    out = cell(nplots,1);

    %% Figure 1: comparison of distance matrices

    hf(1) = figure;
    n    = 1;
    labs = cell(nplots,1);
    for i_g = 1:ngroups
        v.dataFilter.subjectGroup = subject_groups{i_g};

        for i_o = 1:nodor_sets
            thisodorset     = stim_groups{i_o};
            thisodorset_str = stimGroupLabel(thisodorset);
            disp(['Plotting group ''', subject_groups{i_g}, ''' for odors: ', thisodorset_str])

            subplot(ngroups, nodor_sets, n);
            v.dataFilter.stims_allowed = thisodorset;
            out{n} = v.plotDistancesHead('plotType', 'repetitions', 'method', method);

            title(thisodorset_str)
            ylabel(subject_groups{i_g})
            if isempty(out{n}); axis off; end

            labs{n} = [subject_groups{i_g}, ' - ', thisodorset_str];
            n = n + 1;
        end
    end
    set(gcf, 'Position', [1 1 2000 1000])

    %% Figure 2: comparison of distance distributions

    data = [];
    for i = 1:nplots
        if isempty(out{i}); continue; end

        thismat  = 1 - out{i}.distMat3d;
        nslices  = size(thismat, 3);
        thismat  = thismat(repetitions, repetitions, :);
        tri_mask = triu(true(numel(repetitions)), 1);

        out{i}.data = nan(sum(tri_mask, "all"), nslices);
        for j = 1:nslices
            thisslice        = thismat(:,:,j);
            out{i}.data(:,j) = thisslice(tri_mask);
        end

        thisdata = out{i}.data(:);
        data     = [data; i * ones(size(thisdata)), thisdata]; %#ok
    end

    cfg   = v.plotConfig;
    hf(2) = figure;
    boxplot(data(:,2), labs(data(:,1)), 'Orientation', 'horizontal', 'PlotStyle', 'compact', 'Colors', cfg.textcol)
    box off
    set(gca, 'color', cfg.bgcol, 'XColor', cfg.axcol, 'YColor', cfg.axcol);
    set(gcf, 'color', cfg.bgcol);
    xlabel(method, 'Color', cfg.textcol)
    set(gcf, 'Position', [1 1 2000 1000])
    hold off
end

function label = stimGroupLabel(group)
    if iscell(group)
        label = strjoin(group, ', ');
    else
        label = group;
    end
end
