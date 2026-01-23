
function [hf, out, groups, odor_sets] = plotRepetitionDistances(v, method)  
    arguments
        v ExperimentViewer
        method = 'correlation'
    end
    
    % Knobs
    groups = {'naïve', ...
                'trained', ...
                'uncoupled'};
                % 'trained1', ...
                % 'trained2', ...
                % 'trained1 (T-R-S-H-A-ACSF/L)', ...
    % odor_sets = {'all stimuli', ...
    %                 {'Arg','Ala','His','Trp','Ser'}, ...
    %                 {'Leu'}, ...
    %                 'all CS+', ...
    %                 'all CS-', ...
    %                 'all familiar', ...
    %                 'all novel', ...
    %             };
    odor_sets = {{'Arg','Ala','His'}, ...
                    {'Trp','Ser','Leu'}, ...
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
            v.dataFilter.stims_allowed = thisodorset;
            out{n} = v.plotDistancesHead('plotType', 'repetitions', ...
                            'method',method);

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
        nslices = size(thismat,3);
        
        % knobs # TODO : tunable param
        repetitions = 1:5;

        % crop and take out the diagonal and lower triangular
        % matrix
        thismat = thismat(repetitions,repetitions,:);
        idx = triu(true(numel(repetitions)), 1); % only upper triangle idx

        out{i}.data = nan(sum(idx,"all"),nslices);
        for j = 1:nslices
            thisslice = thismat(:,:,j);            
            % return
            out{i}.data(:,j) = thisslice(idx); % column vector
        end
        
        % idx = repmat(idx,1,1,size(thismat,3));
        thisdata = out{i}.data(:);

        % store
        plot_idx = i * ones(size(thisdata));
        data = [data; plot_idx, thisdata];

    end
    
    % plot
    cfg = v.plotConfig;
    hf(2) = figure;
    boxplot(data(:,2),labs(data(:,1)),'Orientation','horizontal','PlotStyle','compact','Colors',cfg.textcol)
    box off
    set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol);
    set(gcf, 'color', cfg.bgcol); 
    xlabel(method,'Color',cfg.textcol)
    set(gcf,'Position',[1 1 2000 1000])
    hold off
end
