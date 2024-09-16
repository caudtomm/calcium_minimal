classdef ActivityTraces
    properties (SetAccess = private)
        PMToff Movie % typically a Snippet
        noLight Movie % typically a Snippet

        ROImap double
        Nrois double
        neuron_IDs double
        background_IDs double
    end
    properties (SetAccess=immutable)
        ntrials double = 0
        t double = []
        N double = 0
        T double = 0
        L double = 0 % single-valued trial duration [frames] (from subject)
        framerate double = 0

        subject_ID char
        source_folder char = pwd
    end
    properties % so called 'derivative properties', because they are derived from everything else
        techBase double = 0
        techNoise double = 0 % single-valued technical noise (variance of PMToff)
        techBaseMap double = [] % temporal avg of technical baseline (expected uniform)
        techNoiseMap double = [] % temporal avg of technical noise (expected uniform)
        
        Fnoise double % fluorescence noise for each ROI [t, roi, trial]
        F0 double % F0 for each ROI, per trial [roi, trial]
        Fpx cell % each {roi} is a [t, px, trials] matrix of raw intensity values
        F double % average F per cell [t, roi, trials] (techBase removed)

        dFoverF double
        
        % currently not implemented #### TODO in separate ActivityTraces = Process(ActivityTraces)
        baseline_periods cell % periods of inactivity for each ROI
        centersurround
        bubble_removal
    end

    methods (Access = private, Static)
        % single-valued number of units
        function N = extractN(subject)
            N = numel(unique(subject.ROImap)) - 1; 
        end
    end

    methods (Access=private)
        % define t [seconds] : 1d array, numel = frames
        function t = extractt(obj,subject)
            fs = subject.framerate;
            t = 0 : 1/fs : obj.L/fs-1/fs;
        end

        % single-valued trial duration [s]
        function T = extractT(obj, subject)
            T = obj.L / subject.framerate;
        end

        function obj = assignROIs(obj,subject)
            neuronmap = subject.ROImap;
            backgroundmap = subject.backgroundROImap;

            bg_ids = unique(backgroundmap);
            obj.background_IDs = bg_ids(bg_ids~=0);

            offset = numel(obj.background_IDs);

            bg_ns = unique(neuronmap);
            bg_ns = bg_ns(bg_ns~=0);
            obj.neuron_IDs = bg_ns + offset;

            neuronmask = ismember(neuronmap,bg_ns);
            neuronmap = neuronmap + neuronmask*offset;
            
            obj.ROImap = backgroundmap + neuronmap;
            obj.Nrois = numel(obj.neuron_IDs) ...
                        + numel(obj.background_IDs);
        end

        % retrieve the baseline intensity value given by technical
        % constraints (ex. PMT noise)
        function [techBase, techBaseMap] = defineTechnicalBaseline(obj)
            techBaseMap = obj.noLight.timeavg;
            techBase = mean(techBaseMap,"all",'omitmissing');
        end

        % retrieve the baseline noise value given by technical
        % constraints (ex. PMT noise variance)
        function [techNoise, techNoiseMap] = defineTechnicalNoise(obj)
            techNoiseMap = var(obj.noLight.stack,[],3);
            techNoise = var(obj.noLight.stack,[],'all');
        end
        
        % F0 for each ROI, per trial [roi, trial]
        function F0 = defineF0(obj)
            % compute (lowest 10th percentile of each ROI's per trial)
            F0 = squeeze(quantile(obj.F,.1,1)); % [roi, trials]
        end
        
        % fluorescence noise for each ROI [t, roi, trial]
        function Fnoise = defineFnoise(obj)
            % initialize output
            Fnoise = nan(obj.L,obj.Nrois,obj.ntrials);

            % compute
            for i_roi = 1:obj.Nrois
                thisroi = obj.Fpx{i_roi}; % [t, px, trials]
                
                % remove technical baseline
                thisroi = thisroi - obj.techBase;

                % compute instantaneouse fluorescence noise as the 
                % variance across ROI px
                noisevals = var(thisroi,[],2,'omitnan'); % [t, 1, trials]
                Fnoise(:,i_roi,:) = noisevals;
            end
        end
        
        % Function to extract raw traces, separately for each pixel, within
        % each ROI. Output Fpx is a cell array when each element is an ROI
        % and contains an [t, px, trials] matrix of raw intensity values
        function Fpx = extractFpx(obj, subject)
            % find movies in the current folder
            movies = cell(obj.ntrials,1);
            for i = 1:obj.ntrials
                movies{i} = find_daughter_file(subject.filelist(i).name,'mat');
            end

            % initialize Fpx output
            Fpx = cell(obj.Nrois,1);

            % find ROI indices (exclude 0 = background)
            rois = unique(obj.ROImap); rois(rois==0) = [];

            % load each movie and extract each ROI's pixel values
            for i_trial = 1:obj.ntrials
                disp(movies{i_trial})
                load(movies{i_trial},'movie')
                for i_roi = 1:obj.Nrois
                    % isolate ROI in the current trial
                    mask = obj.ROImap == rois(i_roi);
                    mask = repmat(mask,1,1,movie.nfr);
                    thisroi = movie.stack(mask);
                    npx = numel(thisroi)/movie.nfr;
                    thisroi = reshape(thisroi,npx,movie.nfr)'; % output [t, px]

                    % scanimage will output a -1 in rare cases of
                    % measurement error. replace those with NaNs
                    thisroi(thisroi==-1) = nan;

                    % initialize the Fpx output for this ROI, according
                    % to the number of pixels it contains
                    if i_trial == 1
                        Fpx{i_roi} = nan(obj.L, npx, obj.ntrials); % obj.L and movie.nfr should be identical!
                    end

                    % store to output cell array
                    Fpx{i_roi}(:,:,i_trial) = thisroi;
                end
            end
        end

        % average F per cell [t, roi, trials] (techBase removed)
        function F = defineF(obj)
            % initialize output
            F = nan(obj.L,obj.Nrois,obj.ntrials);
            
            % fill in each roi's F(t)
            for i = 1:obj.Nrois
                traces = mean(obj.Fpx{i},2,'omitmissing'); % [t, 1, trials]

                % remove technical baseline
                traces = traces - obj.techBase;
                
                % store to output
                F(:,i,:) = traces;
            end
        end

        % dF/F0 per cell [t, roi, trials]
        function dFoverF = definedFoverF(obj)
            % prepare F0 matrix for algebra
            F0rep = shiftdim(obj.F0,-1); % [1,roi,trials]
            F0rep = repmat(F0rep,[obj.L,1,1]);

            % compute
            dFoverF = (obj.F - F0rep) ./ F0rep;

            % get average dFoverF for background ROIs
            background_dFoverF = mean(dFoverF(:,obj.background_IDs,:),2,'omitmissing'); % [t,1,trials]
            background_dFoverF_rep = repmat(background_dFoverF,[1,obj.Nrois,1]);

            % subtract background average from all dFoverF (this is
            % intended to remove global artefacts and fluctuations)
            dFoverF = dFoverF - background_dFoverF_rep;
        end
        
        % SKETCH
        function [centeravg, surroundavg, centermovie, surroundmovie] = ...
                defineCenterSurround(obj)

        end

        % SKETCH
        function traces = desurround(obj)

        end

        % SKETCH
        function traces = removeBubbles(obj)

        end
    end

    methods
        function obj = ActivityTraces(subject) % Constructor
            arguments
                subject Subject
            end

            % define all immutable properties 
            % (determined entirely by the Subject)
            obj.subject_ID = subject.id;
            obj.source_folder = fullfile(subject.locations.subject_datapath,...
                subject.locations.traces_src); 
            obj.framerate = subject.framerate;
            obj.N = obj.extractN(subject);
            obj.L = subject.getNFrames;
            obj.T = obj.extractT(subject);
            obj.t = obj.extractt(subject);
            obj.ntrials = subject.getNTrials;

            obj = assignROIs(obj,subject);
        end

        % check behavior with partial arguments ################
        function obj = loadMovieData(obj,subject,fname_PMToffmeta,fname_noLightmeta)
            arguments
                obj
                subject Subject
                fname_PMToffmeta = ''
                fname_noLightmeta = ''
            end

            % load technical baseline and noise movie Snippets
            thisPMToff = subject.retrieveSnippetfromCoords(fname_PMToffmeta);
            obj = obj.setPMToff(thisPMToff);
            thisnoLight = subject.retrieveSnippetfromCoords(fname_noLightmeta);
            obj = obj.setNoLight(thisnoLight);

            % obj = obj.defineFundamentalProperties(subject);
        end

        function obj = defineFundamentalProperties(obj,subject)
            arguments
                obj
                subject Subject
            end

            % define technical baseline and noise
            [obj.techBase, obj.techBaseMap] = obj.defineTechnicalBaseline();
            [obj.techNoise, obj.techNoiseMap] = obj.defineTechnicalNoise();
            
            % extract pixelwise raw intensity values for each ROI
            obj.Fpx = obj.extractFpx(subject);
        end

        function obj = setDerivativeProperties(obj)
            obj.Fnoise = obj.defineFnoise();
            obj.F = obj.defineF();
            obj.F0 = obj.defineF0();

            obj.dFoverF = obj.definedFoverF();
        end

        % Setter for read-only property 'PMToff' (type Movie)
        function obj = setPMToff(obj, newval)
            arguments
                obj
                newval Movie = Movie() % typically a Snippet
            end
            
            obj.PMToff = newval;
        end

        % Setter for read-only property 'noLight' (type Movie)
        function obj = setNoLight(obj, newval)
            arguments
                obj
                newval Movie = Movie() % typically a Snippet
            end

            obj.noLight = newval;
        end
        
    end

end
