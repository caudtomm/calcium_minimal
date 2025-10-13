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

         function [t_resampled, resampledtrace] = resample(trace,trials,T,L)
            if nargin < 4
                L = median(1+diff(trials,[],2).frame_start); % median trial length [frames]
            end

            s2 = cell(height(trials), 1);
            for i = 1:height(trials)
                s2{i} = trace.raw(trials.frame_start(i):trials.frame_end(i));
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

         function freq = findDominantFrequency(tracemat,fs)
            A = tracemat(:);

            % instantaneous frequency,
            % rolling average over 500 ms
            x = movmean(instfreq(A,fs,'Method','hilbert'),fs/2);

            freq = reshape([x; x(end)],size(tracemat));
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

            % read traces from FiJI output CSVs or summary MATs
            obj.Breathing = obj.extractTraces(fullfiletol(fpath,'rot','Head_vals'));
            obj.Tail = obj.extractTraces(fullfiletol(fpath,'rot','Tail_vals'));

            % highpass the Breathing trace at 0.5 Hz
            obj.Breathing.raw = highpass(obj.Breathing.raw,.5,obj.framerate);

            % clean from scanning background
            MaiTai_freq = 7.66;                                                 % ### TODO: should be argument (specify from Subject)
            obj.Breathing = obj.removeBackground(obj.Breathing,MaiTai_freq);
            obj.Tail = obj.removeBackground(obj.Tail,MaiTai_freq);


            % resample traces to common time axis (LED trials)
            [obj.LED.t_resampled,obj.LED.resampled] = ...
                obj.resample(obj.LED,obj.LED.trials,170);
            [obj.Breathing.t_resampled,obj.Breathing.resampled] = ...
                obj.resample(obj.Breathing,obj.LED.trials,170);
            [obj.Tail.t_resampled,obj.Tail.resampled] = ...
                obj.resample(obj.Tail,obj.LED.trials,170);

            % resample traces to 2p framerate
            [obj.LED.t_resampled2p,obj.LED.resampled2p] = ...
                obj.resample(obj.LED,obj.LED.trials,170,1300);
            [obj.Breathing.t_resampled2p,obj.Breathing.resampled2p] = ...
                obj.resample(obj.Breathing,obj.LED.trials,170,1300);
            [obj.Tail.t_resampled2p,obj.Tail.resampled2p] = ...
                obj.resample(obj.Tail,obj.LED.trials,170,1300);


            % get freq
            obj.Breathing.freq = obj.findDominantFrequency(obj.Breathing.resampled,obj.Breathing.fs);
            obj.Breathing.freq2p = obj.findDominantFrequency(obj.Breathing.resampled2p,MaiTai_freq);
            
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