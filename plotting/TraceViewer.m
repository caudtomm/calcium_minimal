classdef TraceViewer
    properties
        traces ActivityTraces
        idx logical % which units to use
    end

    methods (Static)
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

        function h = plotAvgCurve(obj, t, M)
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

            % plotting
            h = figure; hold on
            b = patch([t,fliplr(t)],[curve1; fliplr(curve2')'],'w','FaceAlpha',.3,'EdgeColor','none');
            b.Annotation.LegendInformation.IconDisplayStyle = 'off';
            plot(t,avgCurve, 'LineWidth', 5, 'Color', 'w')

        end

        function h = plotTrialAvgs(obj, tracename)
            % plots average trace over all cells, for each trial.
            % Staggers trials. Default traces are dFoverF. Assumes input
            % traces : [t, roi, trials]
            arguments
                obj
                tracename char = 'dFoverF'
            end
            
            v = obj.initializeTraces(tracename); % [t, roi, trial]
            v = squeeze(mean(v,2,'omitmissing')); % [t, trial]
            
            h = obj.plotStaggered(obj.traces.t, v);
            xlabel('time [s]')
            ylabel('trial #')
        end

        function h = plotTracesInTime(obj, tracename, trials)
            % plots single-cell traces side-by-side.
            % Staggers cells. Default traces are dFoverF. Assumes input
            % traces : [t, roi, trials]
            arguments
                obj
                tracename char = 'dFoverF'
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

        function h = plotPSTH(obj, tracename, trials, stim_frame, ps_lim)
            arguments
                obj 
                tracename char = 'dFoverF'
                trials double = [1:obj.traces.ntrials]
                stim_frame double = obj.traces.stim_series.frame_onset(1)
                ps_lim double = [-4, 30] % from stim onset [s]
            end

            M = obj.initializeTraces(tracename); % [t, roi, trial]

            fs = obj.traces.framerate;

            interval = floor(ps_lim(1)*fs+stim_frame) : floor(ps_lim(2)*fs+stim_frame); 
            t = interval /fs;

            h = obj.plotAvgCurve(t, M(interval,:,trials));
            xlabel('time from stimulus onset [s]')
            ylabel(tracename)
            axis tight
            
            set(gcf, 'Position', [50 50 400 280]);
            set(gca, 'color', 'none', 'XColor','w', 'YColor','w', 'ZColor','w');
            set(gcf, 'color', 'none'); 
            
            l = legend('PSTH'); legend('boxoff')
            l.TextColor = 'w';
        end

        function h = plotAvgResponse(obj)
            
        end

    end

end