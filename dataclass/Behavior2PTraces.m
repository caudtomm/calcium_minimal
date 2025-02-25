classdef Behavior2PTraces
    properties
        subject_ID char
        framerate double
        LED
        Breathing
        Tail
    end
    
    methods (Static)
        function plotSpectrogram(trace)
            a = trace.raw;
            fs = trace.fs;

            [b,c,d] = pspectrum(a,fs,'spectrogram','FrequencyResolution',1);
            figure; imagesc(d,flipud(c),b)
            yticklabels(flipud(yticklabels));
            title('Spectrogram')
            xlabel('time [s]')
            ylabel('frequency [Hz]')
            colb = colorbar;
            colb.Label.String = 'Power';
        end


        function plotInstRate(trace)
            x = trace.freq;
            t = trace.t;

            figure;plot(t(2:end),x)
            title('Estimated istantaneous rate (rolling avg: 0.5 s)')
            xlabel('time [s]')
            ylabel('frequency [Hz]')
            axis tight
        end

         function plotRawTrace(trace)
            x = trace.raw;
            t = trace.t;

            figure;plot(t,x)
            title('Raw trace')
            xlabel('time [s]')
            ylabel('frequency [Hz]')
            axis tight
         end

         function traceout = extractLEDframestamps(trace)
            % this only makes sense for the odor-system LED
            % initialize output
            traceout = trace;

            % convenience assignments
            x = trace.raw;
            fs = trace.fs;

            %% extract trial periods ('light' refers to the 2p laser, not the LED!)
            
            % get light-on intervals
            th = quantile(x,.001)*2/3; % threshold intensity val (assumes centered data)
            onperiods = convertPeriods(x>th);
            % prune start and end
            if x(1)>th
                onperiods(1,:)=[];
            end
            if x(end)>th
                onperiods(end,:)=[];
            end
            % prune periods under 1 sec (probably artefacts)
            dur = diff(onperiods,[],2);
            onperiods(dur<1*fs,:) = []; 
            % store
            ntrials = size(onperiods,1);
            ITI = median(diff(onperiods(:,1))); % [frames]
            onperiods2p = onperiods; % to keep live until storing

            %% extract LED on periods (...)

            % get light-on intervals
            th = quantile(x,.999)*2/3; % threshold intensity val (assumes centered data)
            onperiods = convertPeriods(x>th);
            
            % prune too long (true odor delivery blips are short!)
            dur = diff(onperiods,[],2);
            onperiods(dur>5,:) = []; 

            % ##################### CONTINUE HERE ######################

            %% store to output
            trials = table();
            trials.frame_start = onperiods2p(:,1);
            trials.frame_end = onperiods2p(:,2);
            
            
            traceout.trials = trials;
            traceout.ntrials = ntrials;
            traceout.ITI = ITI; % [frames]

         end

         function trace = resample(trace,trials,T)
            s2 = cell(height(trials), 1);
            for i = 1:height(trials)
                s2{i} = trace.raw(trials.frame_start(i):trials.frame_end(i));
            end 

            t_common = linspace(0,T,median(1+diff(trials,[],2).frame_start));

            s2_resampled = cell(height(trials),1);
            for i = 1:height(trials)
                t = linspace(0,T,numel(s2{i}));
                s2_resampled{i} = interp1(t,s2{i},t_common,'linear');
            end
            s2_resampled = cell2mat(s2_resampled)';

            trace.t_resampled = t_common;
            trace.resampled = s2_resampled;
         end

         function trace = removeBackground(trace,background_freq)
             s0 = trace.raw;
             
             % s = s0(:);
             % s1 = background_source(:);
             % s_denoised = s - (s1 * (s' * s1) / (s1' * s1));

             n_harmonics = 2;

             s_denoised = s0;
             for i = 1:n_harmonics+1
                noise = bandpass(s0,i*background_freq+[-.3 .3],trace.fs);
                s_denoised = s_denoised-noise;
             end

             trace.raw = s_denoised;

         end

         function freq = findDominantFrequency(trace)
            A = trace.resampled(:);

            % instantaneous frequency,
            % rolling average over 500 ms
            x = movmean(instfreq(A,trace.fs,'Method','hilbert'),trace.fs/2);

            freq = reshape([x; x(end)],size(trace.resampled));
         end
    end

    methods
        function obj = Behavior2PTraces(fpath, framerate)
            arguments
                fpath char = fullfiletol(pwd,'tail_movies')
                framerate double = 50; % Hz
            end

            if isunix
                obj.subject_ID = fliplr(extractBefore(fliplr(pwd),'/'));
            else
                obj.subject_ID = fliplr(extractBefore(fliplr(pwd),'\'));
            end

            obj.framerate = framerate;

            obj.LED = obj.extractTraces(fullfiletol(fpath,'LED_vals'));
            obj.Breathing = obj.extractTraces(fullfiletol(fpath,'Lip_vals'));
            obj.Tail = obj.extractTraces(fullfiletol(fpath,'Tail_vals'));

            % highpass the Breathing trace at 0.5 Hz
            obj.Breathing.raw = highpass(obj.Breathing.raw,.5,obj.framerate);

            % clean from scanning background
            MaiTai_freq = 7.66;
            obj.Breathing = obj.removeBackground(obj.Breathing,MaiTai_freq);
            obj.Tail = obj.removeBackground(obj.Tail,MaiTai_freq);

            % extract framestamps from the LED frame
            obj.LED = obj.extractLEDframestamps(obj.LED);

            % resample traces
            obj.LED = obj.resample(obj.LED,obj.LED.trials,170);
            obj.Breathing = obj.resample(obj.Breathing,obj.LED.trials,170);
            obj.Tail = obj.resample(obj.Tail,obj.LED.trials,170);

            % get freq
            obj.Breathing.freq = obj.findDominantFrequency(obj.Breathing);
        end

        function struct_out = extractTraces(obj,fpath)
            disp(fpath);
            currentDir = pwd;
            cd(fpath)

            % convenience
            fs = obj.framerate;
            
            % get traces
            files = dir('*.csv');
            A = [];
            for ff = 1:numel(files)
                file = files(ff);
                fprintf(['Analysing: ', file.name, '...'])
                
                a = readtable(file.name); a = a.Mean;
                A = [A; a-median(a)];
                
                fprintf(' DONE!\n')
            
            end

            % time axis
            t = [0:1/fs:length(A)/fs-1/fs];

            % store to output
            struct_out.t = t; % time axis
            struct_out.raw = A; % raw intensity trace (median-centered)
            struct_out.fs = fs; % framerate copy

            % return to original directory
            cd(currentDir)
        end

       
    end

end