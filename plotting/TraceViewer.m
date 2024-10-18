classdef TraceViewer
    properties
        traces ActivityTraces
        idx logical % which traces to use
    end

    methods (Static)
        function v = catTraces(M)
            % concatenates Traces:
            % [t , roi , trial] -> [t x trial , roi]
            arguments
                M double
            end

            [T,N,ntrials] = size(M);

            v = nan(T*ntrials,N);
            
            for i = 1:ntrials
                interval = [1:T] + T*(i-1);
                v(interval,:) = M(:,:,i);
            end
        end
    end

    methods (Static, Access=protected)
        function h = plotStaggered(t,M)
            % plots the columns of M staggered in Y

            [nrows,ncols] = size(M);
            
            try
                scaling_fact = mad(M,1,"all");
            catch
                scaling_fact = std(M,[],"all",'omitmissing');
            end
            M = M./scaling_fact;

            % space out data
            spacer_gain = 2; % (2 std)
            spacer = repelem(spacer_gain,ncols);
            spacermat = repmat(cumsum(spacer),[nrows,1]);

            y = M + spacermat;

            % plotting
            h = figure;
            plot(t,y)

            xlim([0 max(t)])
            ylim([0 max(y,[],'all','omitmissing')+spacer_gain])
            yticklabels(yticks./spacer_gain)

            box off


        end
    end

    methods  (Access=protected)
        function v = initializeTraces(obj, tracename)
            arguments
                obj
                tracename = 'dFoverF'
            end
            
            % define which traces to use
            v = eval(strcat('obj.traces.',tracename));
            % select traces of interest only
            v = v(:,obj.idx,:);
        end
    end

    methods
        function obj = TraceViewer(traces, rois_touse)
            arguments
                traces ActivityTraces
                rois_touse double = traces.goodNeuron_IDs
            end
            obj.traces = traces;


            % convert input roi names to use into boolean indices
            allrois = unique(traces.ROImap); allrois = allrois(allrois ~= 0);
            obj.idx = ismember(allrois,rois_touse);
        end

        function plotTrialAvgs(obj, tracename)
            % plots average trace over all cells, for each trial.
            % Staggers trials. Default traces are dFoverF. Assumes input
            % traces : [t, roi, trials]
            arguments
                obj
                tracename = 'dFoverF'
            end
            
            v = obj.initializeTraces(tracename); % [t, roi, trial]
            v = squeeze(mean(v,2,'omitmissing')); % [t, trial]
            
            h = obj.plotStaggered(obj.traces.t, v);
            xlabel('time [s]')
            ylabel('trial #')
        end

        function plotTracesInTime(obj, tracename, trials)
            % plots single-cell traces side-by-side.
            % Staggers cells. Default traces are dFoverF. Assumes input
            % traces : [t, roi, trials]
            arguments
                obj
                tracename = 'dFoverF'
                trials double = [1:obj.traces.ntrials]
            end
            
            M = obj.initializeTraces(tracename); % [t, roi, trial]
            M = M(:,:,trials);
            v = obj.catTraces(M);
            t = linspace(0,size(M,3)*obj.traces.T,size(v,1));

            h = obj.plotStaggered(t, v);
            xlabel('time [s]')
            ylabel('neuron #')
        end

    end

end