classdef Behavior2PTraces
    properties
        subject_ID char
        framerate double
        fpath char
        regions = {'LED', 'Lip', 'Tail'}
        f_range = [] % frame range to analyze [frames]
        ledROI
        lipROI
        tailROI
        LED
        Breathing
        Tail
        headPCA struct
        respiration_PC % index of the PC used for breathing trace (selected manually)
    end

    methods (Access = private)
        function [A,t,fs] = readTraces_csv(obj)
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
            fs = obj.framerate;
            t = [0:1/fs:length(A)/fs-1/fs];
        end

        function [A,t,fs,frange] = readTraces_mat(obj)
            files = dir('*avg_intensities.mat');
            if isempty(files)
                error('No *avg_intensity.mat files found.');
            end
            file = files(1);
            fprintf(['Loading: ', file.name, '...'])
            
            data = load(file.name).avg_intensities;

            frange = obj.f_range;
            if isfield(data, 'allchecked')
                A = data.allchecked;
                frange = 1:length(A);
            elseif isfield(data, 'all')
                A = data.all;
                if ~isempty(frange)
                    A = A(frange);
                end
            else
                error('The loaded file does not contain the field "all".');
            end

            A = A-median(A,'omitmissing');                
            A = fillmissing(A,'linear');
            
            fprintf(' DONE!\n')

            % time axis
            fs = obj.framerate;
            t = [0:1/fs:length(A)/fs-1/fs];
        end
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

            %% inference of camera frame rate from measured ITI vs expected ITI
            expected_ITI_s = 180; % [s] (assumes 3 min ITI)
            inferred_framerate = ITI/expected_ITI_s; % [Hz]
            fprintf('Inferred camera framerate: %.2f Hz\n',inferred_framerate);

            %% extract LED-on periods (...)

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
            traceout.fs = inferred_framerate; % update framerate

         end

         function [t_resampled, resampledtrace] = resample(tracemat,trials,T,L)
            if nargin < 4
                L = median(1+diff(trials,[],2).frame_start); % median trial length [frames]
            end

            rawtrace = tracemat(:);

            s2 = cell(height(trials), 1);
            for i = 1:height(trials)
                s2{i} = rawtrace(trials.frame_start(i):trials.frame_end(i));
            end 

            t_common = linspace(0,T,L);

            s2_resampled = cell(height(trials),1);
            for i = 1:height(trials)
                t = linspace(0,T,numel(s2{i}));
                s2_resampled{i} = interp1(t,s2{i},t_common,'linear');
            end
            s2_resampled = cell2mat(s2_resampled)';

            t_resampled = t_common;
            resampledtrace = s2_resampled;
         end

         function trace = removeBackground(trace,fs,background_freq, method)
             s0 = trace(:);

             n_harmonics = 2;

            if nargin < 4
                method = 'notch';
            end

            switch method
                case 'notch'
                    s_denoised = notch_harmonics(s0, fs, background_freq, 0.25, 0.95);
                case 'bandpass'
                    s_denoised = bandpass_harmonics(s0, fs, background_freq, .25, n_harmonics);
                otherwise
                    error('Unknown method for background removal.');
            end

            trace = s_denoised;


            function s_denoised = bandpass_harmonics(s0, fs, background_freq, bw, n_harmonics)
                % iterative bandpass filtering
                % DEPRECATED: use notch instead
                % (can generate artifacts from sequential bandpass filtering)
                % (doesn't care about nyquist - can be bad for low fs)
                s_denoised = s0;
                for i = 1:n_harmonics+1
                    noise = bandpass(s0,i*background_freq+[-bw bw],fs); % zero-phase
                    s_denoised = s_denoised-noise;
                end
            end

            function y = notch_harmonics(x, fs, f0, bw, steep)
                if nargin<5, steep = 0.95; end
                y = x;
                K = floor((fs/2 - 1e-6)/f0);       % harmonics below Nyquist
                for k = 1:K
                    f = k*f0;
                    y = bandstop(y,[max(0,f-bw) min(fs/2-1e-6,f+bw)],fs,'Steepness',steep); % zero-phase
                end
            end
         end

         function freq = findDominantFrequency(tracemat,fs, method) % # DEPRECATED
            A = tracemat(:);

            if nargin < 3
                method = 'hilbert';
            end

            % compute the instantaneous frequency

            switch method
                case 'hilbert'
                    % Hilbert transform method
                    x = movmean(instfreq(A,fs,'Method','hilbert'),fs/2); % rolling average over 500 ms
                    freq = reshape([x; x(end)],size(tracemat));
                otherwise
                    error('Unknown method for instantaneous frequency estimation.');
            end
            
         end

         function visualizeBreathing(traceStruct)
            % VISUALIZEBREATHING  Quick visual QC for breathing event detection.
            % Input:
            %   traceStruct.resampled   - raw respiration trace
            %   traceStruct.eventsFP    - event timestamps (s)
            %   traceStruct.rateFP      - instantaneous rate (Hz)
            %   traceStruct.fs          - sampling frequency (Hz)

            raw = traceStruct.resampled(:);
            fs = traceStruct.fs;
            t = (0:numel(raw)-1)/fs;

            figure('Name','Breathing QC','Color','w','Position',[100 100 1000 500])

            subplot(2,1,1)
            hold on
            plot(t,raw,'Color',[0.7 0.7 0.7])
            if isfield(traceStruct,'eventsFP') && ~isempty(traceStruct.eventsFP)
                %xline(find(traceStruct.eventsFP(:))./fs)
            end
            xlabel('Time (s)')
            ylabel('Amplitude')
            title('Raw and Processed Traces with Detected EvSents')
            legend({'Raw','Processed','Events'},'Location','best')
            axis tight
            ax(1) = gca;

            subplot(2,1,2)
            if isfield(traceStruct,'rateFP') && ~isempty(traceStruct.rateFP)
                plot(t, traceStruct.rateFP(:),'k','LineWidth',1.2)
                xlabel('Time (s)')
                ylabel('Rate (Hz)')
                title('Instantaneous Breathing Rate')
                axis tight
            else
                text(0.5,0.5,'No rate data available','HorizontalAlignment','center')
                axis off
            end
            ax(2) = gca;

            linkaxes(ax,'x')
            sgtitle('Breathing Signal Overview')

            disp('Press ENTER to continue...')
            pause

            close(gcf)
        end


         function [events, ipi, rate] = findBreathingEvents(tracemat,fs,method,minD)
            if nargin < 3
                method = 'findpeaks';
            end
            trace = tracemat(:);
            if nargin < 4
                minD = round(0.1*fs); % minimum interpeak distance is statically set to 100 ms (10 Hz max breathing rate)
            end

            switch method
                case 'findpeaks'
                    [~,locs] = findpeaks(trace, ...
                        'MinPeakProminence',0.8, ...
                        'MinPeakWidth',round(0.03*fs), ...
                        'MinPeakDistance',minD); % minimum interpeak distance is statically set to 100 ms (10 Hz max breathing rate)
                    events = locs;
                case 'RF'
                    % use Rainer's favorite function
                    events = RF_peakdetect(trace);
                        
                    events(diff([0;events])<minD) = []; % enforce minimum interpeak distance
                otherwise
                    error('Unknown method for breathing event detection.');
            end

            % ensure events are a sorted column vector and within bounds
            events = unique(events(:));
            events(events < 1) = [];
            events(events > numel(trace)) = [];

            n = numel(trace);
            rate = NaN(n,1);

            if numel(events) < 2
                ipi = [];
                return
            end

            % compute inter-peak intervals (s) and instantaneous interval rates (Hz)
            tpk = events ./ fs;
            ipi = diff(tpk);               % seconds
            r_intervals = 1 ./ ipi;        % Hz

            % assign interval rate to samples between consecutive peaks
            for k = 1:numel(r_intervals)
                istart = events(k);
                iend = events(k+1)-1;
                rate(istart:iend) = r_intervals(k);
            end
            % assign last interval rate from last peak to end
            rate(events(end):end) = r_intervals(end);

            % interpolate missing samples over time and smooth (~1 s gaussian)
            t_samples = (0:n-1)' ./ fs;
            if any(~isnan(rate))
                rate = fillmissing(rate,'linear','SamplePoints',t_samples);
                rate = fillmissing(rate,'nearest'); % fill start/end
                rate = smoothdata(rate,'gaussian',round(1*fs));
            end

            % convert to logical array with length = numel(tracemat)
            tmp = false(n,1);
            tmp(events) = true;
            events = tmp;
            % to get indices, use find(events)

            % fold back to original shape
            events = reshape(events,size(tracemat));
            rate = reshape(rate,size(tracemat));

         end
    end

    methods
        function obj = Behavior2PTraces(fpath, ledROI,lipROI,tailROI, varargin)
            arguments
                fpath char = fullfiletol(pwd,'tail_movies') 
                ledROI = []
                lipROI = []
                tailROI = []
            end
            arguments (Repeating)
                varargin
            end

            if isunix
                obj.subject_ID = fliplr(extractBefore(fliplr(pwd),'/'));
            else
                obj.subject_ID = fliplr(extractBefore(fliplr(pwd),'\'));
            end

            % save region coordinates
            obj.ledROI = ledROI;
            obj.lipROI = lipROI;
            obj.tailROI = tailROI;

            % extract stated framerate from video files
            obj.fpath = fpath;
            videos = dir(fullfiletol(fpath,'*.avi'));
            obj.framerate = VideoReader(fullfiletol(fpath,videos(1).name)).FrameRate;

            % generate region crops and save to folders
            %obj.cropMovies

            % extract framestamps from the LED frame
            obj.LED = obj.extractTraces(fullfiletol(fpath,'rot','LED_vals'));
            obj.LED = obj.extractLEDframestamps(obj.LED);
            obj.f_range = obj.LED.f_range; % update frame range to analyze
            obj.framerate = obj.LED.fs; % update framerate

            % read PCA of head motion (for the moment, this is unused here)
            fileIn = fullfiletol(fpath,'rot','Head_vals','pca_results.mat');
            if isfile(fileIn)
                obj.headPCA = load(fileIn);
            else
                warning('PCA results file for head motion not found.');
                obj.headPCA = struct();
            end
            if exist('headPC_fpath','var') && isfile(headPC_fpath)
                obj.respiration_PC = load(headPC_fpath).respiration_PC;
            end

            % read traces from FiJI output CSVs or summary MATs
            obj.Breathing = obj.extractTraces(fullfiletol(fpath,'rot','Head_vals'));
            [obj.Breathing.raw, obj.Breathing.t] = obj.mergeHeadPCtrace; % use PC trace if available
            obj.Tail = obj.extractTraces(fullfiletol(fpath,'rot','Tail_vals'));

            % resample traces to common time axis (LED trials)
            trial_length = 170; % sec
            [obj.LED.t_resampled,obj.LED.resampled] = ...
                obj.resample(obj.LED.raw,obj.LED.trials,170);
            [obj.Breathing.t_resampled,obj.Breathing.resampled] = ...
                obj.resample(obj.Breathing.raw,obj.LED.trials,170);
            [obj.Tail.t_resampled,obj.Tail.resampled] = ...
                obj.resample(obj.Tail.raw,obj.LED.trials,170);

            % temporarily linearize
            breath_dims = size(obj.Breathing.resampled);
            obj.Breathing.resampled = obj.Breathing.resampled(:);
            tail_dims = size(obj.Tail.resampled);
            obj.Tail.resampled = obj.Tail.resampled(:);

            % clean from scanning background
            scanimage_framenum = 1300;
            MaiTai_freq = scanimage_framenum/trial_length;   % ### TODO: should be argument/-s (specify from Subject)
            processed = obj.removeBackground(obj.Tail.resampled,obj.Tail.fs,MaiTai_freq);
            obj.Tail.resampled = reshape(processed,tail_dims); % overwrite resampled with processed for consistency
            % breathing trace will be cleaned later inside processBreathing

            % detect breathing events
            processed = obj.processBreathing(obj.Breathing,MaiTai_freq); 
            processed = reshape(processed,breath_dims);
            obj.Breathing.resampled = processed; % overwrite resampled with processed for consistency
            [obj.Breathing.eventsRF,obj.Breathing.ipiRF,obj.Breathing.rateRF] = ...
                obj.findBreathingEvents(processed,obj.Breathing.fs,'RF');
            [obj.Breathing.eventsFP,obj.Breathing.ipiFP,obj.Breathing.rateFP] = ...
                obj.findBreathingEvents(processed,obj.Breathing.fs,'findpeaks');

            % resample traces to 2p framerate
            [obj.LED.t_resampled2p,obj.LED.resampled2p] = ...
                obj.resample(obj.LED.raw,obj.LED.trials,170,1300);
            resampled_trials = obj.defResampledTrials(size(obj.Breathing.resampled));
            [obj.Breathing.t_resampled2p,obj.Breathing.resampled2p] = ...
                obj.resample(obj.Breathing.resampled(:),resampled_trials,170,1300);
            [obj.Tail.t_resampled2p,obj.Tail.resampled2p] = ...
                obj.resample(obj.Tail.raw,obj.LED.trials,170,1300);

            % get freq
            % obj.Breathing.freq = obj.findDominantFrequency( ...
            %     obj.Breathing.resampled,obj.Breathing.fs, 'hilbert');
            % obj.Breathing.freq2p = obj.findDominantFrequency( ...
            %     obj.Breathing.resampled2p,MaiTai_freq, 'hilbert');

            % get freq (RF)
            tmp = reshape(obj.Breathing.rateRF,[],35);
            tx = linspace(0,170,height(tmp));
            [tmp,ty] = resample(tmp,tx,MaiTai_freq);
            obj.Breathing.rate2p_RF = tmp(2:end,:);
            obj.Breathing.t_rate2p = ty(2:end);

            % get freq (FP)
            tmp = reshape(obj.Breathing.rateFP,[],35);
            tx = linspace(0,170,height(tmp));
            tmp = resample(tmp,tx,MaiTai_freq);
            obj.Breathing.rate2p_FP = tmp(2:end,:);


        end

        function trials = defResampledTrials(obj,sz)
            L = sz(1); ntrials = sz(2);
            frame_start = 1 + ([1:ntrials]'-1).*L;
            frame_end = frame_start + L-1;
            trials = table(frame_start,frame_end);
        end

        function [rawtrace, t] = mergeHeadPCtrace(obj)
            rawtrace = obj.Breathing.raw; % default
            t = obj.Breathing.t; % default
            
            if ~isempty(obj.respiration_PC) && isfield(obj.headPCA,'V')
                % send warnings, but merge anyway
                if obj.respiration_PC < 1 || obj.respiration_PC > size(obj.headPCA.V,2)
                    warning('respiration_PC index (%d) is out of range for headPCA.V (columns: %d).', ...
                        obj.respiration_PC, size(obj.headPCA.V,2));
                else
                    lenV = size(obj.headPCA.V,1);
                    lenB = numel(obj.Breathing.raw);
                    if lenV ~= lenB
                        warning('Length mismatch between headPCA.V (rows=%d) and Breathing.raw (len=%d). PCA vector will be used; ensure traces are aligned or resampled.', ...
                            lenV, lenB);
                    end

                    t = (0:1/obj.framerate:(lenV-1)/obj.framerate);
                end

                rawtrace = obj.headPCA.V(:,obj.respiration_PC);
                
            end
        end

        function processed = processBreathing(obj,trace,MaiTai_freq)
            fs = trace.fs;
            x1 = trace.resampled(:);
            
            % preprocess breathing trace
            processed = highpass(x1,.2,fs); % highpass the Breathing trace at 0.2 Hz (to remove slow drift)
            processed = obj.removeBackground(processed,fs,MaiTai_freq); % remove 2p scanning background up to 2nd harmonic
            
            % bandpass filter to respiration band
            expected_breathing_freq = [.5 12]; % [Hz] % generous band limits
            [b2,a2] = butter(4, expected_breathing_freq/(fs/2));      % respiration band
            processed = filtfilt(b2,a2,processed); % zero-phase

            % detrend and z-score within a sliding window (robust to drifts)
            processed = (processed - movmedian(processed,round(5*fs))) ./ max(eps, movmad(processed,round(5*fs),1));

            % remove filter edge effects
            n = numel(processed);
            m = round(3*fs);                     % discard 3 s from each end
            mask = true(n,1);
            mask(1:m) = false; mask(end-m+1:end) = false;
            processed(~mask) = NaN;

        end

        function struct_out = extractTraces(obj,fpath)
            disp(fpath);
            currentDir = pwd;
            cd(fpath)
            
            % get traces
            %[A,t,fs] = obj.readTraces_csv();
            [A,t,fs,frange] = obj.readTraces_mat();

            % store to output
            struct_out.t = t; % time axis
            struct_out.raw = A; % raw intensity trace (median-centered)
            struct_out.fs = fs; % framerate copy
            struct_out.f_range = frange; % frame range of the trace (empty if not available - unchecked trace)

            % return to original directory
            cd(currentDir)
        end

        function cropMovies(obj)
            disp(obj.fpath);
            currentDir = pwd;
            cd(obj.fpath)
            
            cropnames = obj.regions;
            cropcoords = {obj.ledROI,obj.lipROI,obj.tailROI};

            % get traces
            files = dir('*.avi');
            for ff = 1:numel(files)
                file = files(ff);
                fprintf(['Cropping: ', file.name, '...'])

                for i_crop = 1:numel(cropnames)
                    fout = fullfile(cropnames{i_crop},file.name);
                    save_avi_crop(file.name, fout, cropcoords{i_crop})
                end
                
                fprintf(' DONE!\n')
            
            end

            % return to original directory
            cd(currentDir)
        end
     
    end

end