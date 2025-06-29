classdef TraceViewer
    properties
        traces ActivityTraces
        idx logical % which units to use
    end

    methods (Static)
        function [events, intervals, t] = getPeriEventData(M, event_startfr, ps_lim, fs)
            [L, nrois, nreplicas] = size(M);
            t = ps_lim(1) : 1/fs : ps_lim(2)-1/fs;

            % define peri-event frame intervals
            interval_len = floor(diff(ps_lim)*fs);
            nevents = numel(event_startfr); % number of unique events! amended later, to take replicas into account
            intervals = nan(nevents,interval_len); % [events, frames]
            for i = 1:nevents
                intervals(i,:) = floor(ps_lim(1)*fs+event_startfr(i)) + ...
                    [0 : interval_len-1];
            end
            intervals(intervals>L | intervals<1) = nan;

            % % multiply intervals by number of replicas
            % nevents = nevents * nreplicas; % total number of events
            % ntervals = repmat(intervals,[nreplicas,1]); % [events, frames]

            % isolate peri-event data
            events = nan(interval_len,nrois,nevents*nreplicas); % [t, roi, event]
            for i = 1:nevents
                isnanframe = isnan(intervals(i,:));
                thisevent = M(intervals(i,~isnanframe),:,:);
                
                % actual indexing (taking accont of replicas)
                idx = (i-1)*nreplicas + (1:nreplicas);

                events(~isnanframe,:,idx) = thisevent;
            end
        end

        function v = catTraces(M,dim)
            % concatenates Traces:
            % dim = 1;
            % [t , roi , trial] -> [t x trial , roi]
            % dim = 2;
            % [t , roi , trial] -> [t , roi x trial]
            arguments
                M double
                dim double
            end

            [T,N,ntrials] = size(M);

            switch dim
                case 1
                    v = nan(T*ntrials,N);
                    
                    for i = 1:ntrials
                        interval = [1:T] + T*(i-1);
                        v(interval,:) = M(:,:,i);
                    end
                case 2
                    v = nan(T,ntrials*N);
                    
                    for i = 1:ntrials
                        interval = [1:N] + N*(i-1);
                        v(:,interval) = M(:,:,i);
                    end
            end
        end

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
            spacer = repelem(spacer_gain,ncols) .* std(M,'omitmissing');
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
                tracename = 'dFoverF_good'
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
                rois_touse double = 1:traces.N
            end
            obj.traces = traces;


            % convert input roi names to use into boolean indices
            allrois = unique(traces.ROImap); allrois = allrois(allrois ~= 0);
            obj.idx = ismember(allrois,rois_touse);
        end

        function [h, M, curves] = plotAvgCurve(obj, t, M)
            arguments
                obj 
                t double
                M double
            end

            M = obj.catTraces(M,2);
            
            avgCurve = mean(M,2,'omitmissing');
            
            shade = std(M,[],2,'omitmissing');
            curve1 = avgCurve + shade;
            curve2 = avgCurve - shade;

            % return
            curves = [avgCurve,curve1,curve2];

            % plotting
            h = figure; hold on
            b = patch([t,fliplr(t)],[curve1; fliplr(curve2')'],'w','FaceAlpha',.3,'EdgeColor','none');
            b.Annotation.LegendInformation.IconDisplayStyle = 'off';
            plot(t,avgCurve, 'LineWidth', 5, 'Color', 'w')

        end

        function [h,v] = plotTrialAvgs(obj, tracename)
            % plots average trace over all cells, for each trial.
            % Staggers trials. Default traces are dFoverF. Assumes input
            % traces : [t, roi, trials]
            arguments
                obj
                tracename char = 'dFoverF_good'
            end
            
            v = obj.initializeTraces(tracename); % [t, roi, trial]
            v = squeeze(mean(v,2,'omitmissing')); % [t, trial]
            
            h = obj.plotStaggered(obj.traces.t, v);
            xlabel('time [s]')
            ylabel('trial #')
        end

        function [h,v] = plotTracesInTime(obj, tracename, trials)
            % plots single-cell traces side-by-side.
            % Staggers cells. Default traces are dFoverF. Assumes input
            % traces : [t, roi, trials]
            arguments
                obj
                tracename char = 'dFoverF_good'
                trials double = [1:obj.traces.ntrials]
            end
            
            M = obj.initializeTraces(tracename); % [t, roi, trial]
            M = M(:,:,trials);
            v = obj.catTraces(M,1);
            t = linspace(0,size(M,3)*obj.traces.T,size(v,1));

            h = obj.plotStaggered(t, v);
            xlabel('time [s]')
            ylabel('neuron #')
        end

        function [h, events, plot_data] = plotPSTH(obj, tracename, ps_lim, event_frames, do_compare_trials, event_order)
            arguments
                obj 
                tracename char = 'dFoverF_good'
                ps_lim double = [-4, 30] % from event frame [s]
                event_frames double = [obj.traces.stim_series.frame_onset(1)]
                do_compare_trials logical = true
                event_order double = [];
            end

            M = obj.initializeTraces(tracename); % [t, roi, trial]
            
            if ~do_compare_trials
                % use 2D-formatted data matrix (event timestamps are independent of trial boundaries)
                M = ActivityTraces.format(M,obj.traces.ntrials);
            end

            % From now on, the workflow is the same for both cases. Dim 3 of M is used to "copy" event intervals.
            % Normally, this is useful to compare the same interval across trials (ex. stimulus response). 

            % crop out event data. events : [t, roi, event]
            [events, ~, t] = obj.getPeriEventData(M, event_frames, ps_lim, obj.traces.framerate);

            % sort as specified
            if ~isempty(event_order)
                events = events(:,:,event_order);
            end

            % plot
            [h, plot_data] = obj.plotAvgCurve(t, events);
            xlabel('time from stimulus onset [s]')
            ylabel(tracename)
            axis tight
            
            set(gcf, 'Position', [50 50 400 280]);
            set(gca, 'color', 'none', 'XColor','w', 'YColor','w', 'ZColor','w');
            set(gcf, 'color', 'none'); 
            
            l = legend('PSTH'); legend('boxoff')
            l.TextColor = 'w';
        end

        function [h, events] = plotbehaviorPSTH(obj, tracename, trials, event_frames, ps_lim)
            arguments
                obj 
                tracename char = 'Tail'
                trials double = [1:obj.traces.ntrials]
                event_frames double = [obj.traces.stim_series.trialnum, obj.traces.stim_series.frame_onset]
                ps_lim double = [-4, 30] % from event frame [s]
            end

            M = obj.traces.behavior2p.(tracename).raw;
            

            nrois = size(M,2);
            fs = obj.traces.framerate;
            t = ps_lim(1) : 1/fs : ps_lim(2)-1/fs;

            % define intervals
            interval_len = floor(diff(ps_lim)*fs);
            event_frames = event_frames(ismember(event_frames(:,1),trials),:); % only use events from specified trials
            nevents = size(event_frames,1);
            intervals = nan(nevents,interval_len);
            for i = 1:nevents
                intervals(i,:) = floor(ps_lim(1)*fs+event_frames(i,2)) + ...
                    [0 : interval_len-1];
            end
            intervals(intervals>obj.traces.L | intervals<1) = nan;

            % isolate peri-event data
            events = nan(interval_len,nrois,nevents); % [t, roi, event]
            for i = 1:nevents
                isnanframe = isnan(intervals(i,:));
                thisevent = M(intervals(i,~isnanframe),:,event_frames(i,1));
                
                events(~isnanframe,:,i) = thisevent;
            end

            % plot
            h = obj.plotAvgCurve(t, events);
            xlabel('time from stimulus onset [s]')
            ylabel(tracename)
            axis tight
            
            set(gcf, 'Position', [50 50 400 280]);
            set(gca, 'color', 'none', 'XColor','w', 'YColor','w', 'ZColor','w');
            set(gcf, 'color', 'none'); 
            
            l = legend('PSTH'); legend('boxoff')
            l.TextColor = 'w';
        end

    end

end