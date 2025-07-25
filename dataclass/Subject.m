classdef Subject
    properties
        locations Locations = Locations
        id
        name
        log cell = {}
        group
        traces ActivityTraces
        odor_delay double {mustBeNonnegative}
        notes string = "" 

        reference_img double
        reference_img_meta
        singletrial_meta cell
        anatomy_imgs double
        localcorr_imgs double

        filelist
        trialnum
        min_trialnum
        stim_series
        ch2stimtype_map

        % fspecs
        imaging_date
        owner char
        tg_line char
        withinday_id char
        region char
        protocol char
        method char
        extra_info char

        ROImap double = []
        backgroundROImap double = []
        framerate
        scanimage_metadata
        badperiods
        anat_regions
        
        process_log cell = {}
        currentstate
    end

    properties (SetAccess = protected)
        badtrials double = []
    end

    methods (Access = protected)
        % Method to set list of files
        function [files,trialnum,min_trialnum] = setFileList(obj)
            loc = obj.locations;
            relative_pathin = fullfiletol(loc.subject_ID, loc.orig_trials);

            fprintf(strcat('\nSource folder: ', strrep(relative_pathin, '\','\\'),'.\n\n'))

            files = dir(fullfiletol(loc.general_datapath, relative_pathin, ...
                strcat('*.',loc.datafile_ext)));

            % order by trial number
            n = [];
            for i_f = 1:numel(files)
                file = files(i_f);
                fspecs = getFileNameSpecs(file.name);
                n = [n;fspecs.trial_num];
            end
            min_trialnum = min(n); n=n-min_trialnum+1;
            [~,idx] = sort(n); trialnum = n(idx);
            files = files(idx);
        end
        
        % Method to read series of stimuli (output file of LabVIEW odor 
        % system control script)
        function [stim_series, ch2stimtype_map] = setStimSeries(obj)
            loc = obj.locations;
            fname = strcat(loc.subject_ID,'_',num2str(obj.min_trialnum, '%05.f'));
            FileIn = fullfiletol(loc.subject_datapath,loc.orig_trials,fname);
            [stim_series, ch2stimtype, ch2stimtype_map] = read_stim_series(FileIn);
            stim_series = horzcat(array2table(stim_series, ...
                'VariableNames',{'trialnum' 'odor_channel' 'frame_onset','frame_offset'}), ...
                cell2table(ch2stimtype,'VariableNames',{'stimulus'}));
        end

        % Method to read imaging metadata from scanimage-a
        function [framerate, scanimage_metadata] = setScanimageMetadata(obj) % FROM THE FIRST FILE
            FileIn = fullfiletol(obj.filelist(1).folder,obj.filelist(1).name);
            [A,result,framerate,zstep,zoom,motorpositions,scalingfactors] = ...
                read_metadata_function(FileIn);
            movie = double(loadTiffStack(FileIn));
            [height,width,numberframes] = size(movie);
            
            % package metadata into structure
            scanimage_metadata.A = A;
            scanimage_metadata.result = result;
            scanimage_metadata.framerate = framerate;
            scanimage_metadata.zstep = zstep;
            scanimage_metadata.zoom = zoom;
            scanimage_metadata.motorpositions = motorpositions;
            scanimage_metadata.scalingfactors = scalingfactors;
            scanimage_metadata.height = height;
            scanimage_metadata.width = width;
            scanimage_metadata.numberframes = numberframes;
        end

        % Method to define stimulus-triggered timestamps
        function [stim_on,stim_off,f0_window,response_window] = set_stimtrig_timestamps(obj,mode,i_file)
            fs = obj.framerate;
            % obj.odor_delay is in seconds
            
            if ~exist('i_file','var'); i_file = 1; end
            stim_on = obj.stim_series.frame_onset(i_file) + obj.odor_delay*fs;
            stim_off = obj.stim_series.frame_offset(i_file) + obj.odor_delay*fs;
            f0_window = [1, stim_on-10];
            response_window = floor([stim_on,stim_off] + obj.odor_delay*fs);

            switch mode
                case 'frames'
                    return
                case 'seconds'
                    stim_on = max([0, stim_on/fs]);
                    stim_off = max([0, stim_off/fs]);
                    f0_window = [max([0, f0_window(1)/fs]), ...
                        min([obj.scanimage_metadata.numberframes, f0_window(2)/fs])];
                    response_window = [max([0, response_window(1)/fs]), ...
                        min([obj.scanimage_metadata.numberframes, response_window(2)/fs])];
                otherwise
                    error('Specified stimulus-triggered timestamps extraction mode is unknown.')
            end
        end

        
        
        % Method to load a pre-existing annotation for ROI selection
        function p_ann = load_partial_annotation(obj)
            obj.reload;
            complete_files = obj.ROIcheckfiles.is_complete;
            % if none of the subject files are completely ROId, initialize
            % as an empty vector
            p_ann = [];
            if isempty(complete_files); return; end
            
            % if some files have been ROId completely, ROIs should be copied from
            % the latest completely ROId file
            [~, latest_complete_idx] = max(arrayfun(@(x) getFileNameSpecs(x{:}).trial_num, complete_files));
            fname = getFileNameSpecs(complete_files{latest_complete_idx}).fname;
            FileIn = fullfiletol(obj.locations.subject_datapath,obj.locations.defRois,[fname,'_defROIs.mat']);
            plane = robust_io('load',FileIn,'plane').plane;
            p_ann = plane{1}.ROI_map;
        end
        
    end

    methods
        % Constructor
        function obj = Subject(name,locations,group)
            obj.locations = locations;
            obj.id = locations.subject_ID;
            obj.name = name;
            obj.group = group;
            [obj.filelist, obj.trialnum, obj.min_trialnum] = obj.setFileList();
            [obj.stim_series, obj.ch2stimtype_map] = obj.setStimSeries();
            [obj.framerate, obj.scanimage_metadata] = obj.setScanimageMetadata();

            try
                [obj.reference_img, obj.reference_img_meta] = ...
                    obj.loadReferenceImg(obj.locations.references.raw);
            catch
                warning('Reference folder or image not found')
            end

            try
                obj.badperiods = obj.loadBadperiods();
            catch
                warning('Bad periods not found')
            end

            obj.Tiff2Movie(obj.locations.orig_trials);

            % load and save singletrial metadata separately. retain info
            % only for trial 1
            singletrial_meta = obj.loadSingletrialMeta();
            s.singletrial_meta = singletrial_meta;
            robust_io('save',fullfiletol(obj.locations.subject_datapath,'singletrial_meta.mat'),s,'-v7.3');
            obj.singletrial_meta = cell(size(singletrial_meta));
            obj.singletrial_meta(1) = singletrial_meta(1);
            clear singletrial_meta
            
            obj.anatomy_imgs = obj.retrieve_trial_anatomies();
            obj.localcorr_imgs = obj.retrieve_localcorr_maps();
            obj = obj.setImagingInfo;


            obj = update_currentstate(obj, 'newly constructed');
            % obj.ROIcheckfiles = struct('to_keep',[],'to_discard',[],'is_complete',[]);

            obj.save2mat();
        end
        
        %% Setter methods

        % Method to load and set the calcium imaging reference image (1 for
        % each plane, concatenated in Z)
        function [img,img_meta] = loadReferenceImg(obj,subdirectory)
            disp(['Loading reference image from: ', subdirectory]);

            [img,img_meta] = deal([]);
            loc = obj.locations;
            PathIn = fullfiletol(loc.subject_datapath, subdirectory);
            
            % load reference image
            fname = '*.mat';
            files = dir(fullfiletol(PathIn,fname));
            if ~isempty(files)
                switch numel(files)
                    case 1
                    otherwise
                        warning('more than one reference image available. selected: %s', ...
                            files(1).name)
                end
                snip = robust_io('load',fullfiletol(PathIn,files(1).name)).movie;
                img = snip.timeavg;
                img_meta = snip.metadata;
            else
                fname = strcat('*.',loc.datafile_ext);
                files = dir(fullfiletol(PathIn,fname));
                switch numel(files)
                    case 0
                        warn_notfound('raw reference image'); return
                    case 1
                    otherwise
                        warning('more than one reference image available. selected: %s', ...
                            files(1).name)
                end
                img = double(loadTiffStack(char(fullfiletol(PathIn,files(1).name))));
                
                % load metadata
                fname = 'metadata.json';
                img_meta = readJson(fullfiletol(PathIn,fname));
            end
        end

        function obj = setBadtrials(obj,value)
            arguments
                obj 
                value double
            end

            obj.badtrials = value;

            if isempty(obj.traces); return; end

            obj.traces = obj.traces.setBadtrials(value);
        end

        % Method to set the imaging info
        function obj = setImagingInfo(obj)
            str = obj.singletrial_meta{1}.path.date;
            mydate = datestr(datenum(str ,'yyyymmdd'),'dd/mm/yyyy');
            obj.imaging_date = mydate;

            fspecs = getFileNameSpecs(obj.filelist(1).name);
            obj.owner = fspecs.owner;
            obj.tg_line = fspecs.subject_line;
            obj.withinday_id = fspecs.subject;
            obj.region = fspecs.region;
            obj.protocol = fspecs.stim_type;
            obj.method = fspecs.method;
            obj.extra_info = fspecs.extra;
        end
        
        % setter for secundary properties requiring user input
        function obj = setManually(obj)
            obj.anat_regions = selectAnatRegions(obj,false);
        end

        % Method to load all raw trials and return a cell
        % array with all their metadata
        function result = loadSingletrialMeta(obj,sourcedir)
            loc = obj.locations;
            if ~exist('sourcedir','var') || isempty(sourcedir)
                sourcedir = fullfiletol(loc.subject_datapath,loc.orig_trials);
            end
            files = obj.filelist;

            result = {};

            parfor i_f = 1:numel(files)
                thisfilename = find_daughter_file(fullfiletol(sourcedir,files(i_f).name),'mat');
                if ~isempty(thisfilename)
                    % thismovie = load(thisfilename,'movie').movie;
                    thismovie = robust_io('load',thisfilename,'movie').movie;
                else
                    thisfilename = fullfiletol(sourcedir,files(i_f).name);
                    thismovie = Movie(thisfilename);
                end

                thismeta = objectToPropsStruct(thismovie, {'stack'}); % stack is excluded because too heavy
                result{i_f} = thismeta;
            end
        end
        
        % Method to load BADPERIOD coordinates
        function badperiods = loadBadperiods(obj)
            loc = obj.locations;
            FileIn = fullfiletol(loc.subject_datapath,loc.orig_trials,[loc.subject_ID,'_badperiods.csv']);
            badperiods.data = [];
            if exist(FileIn,'file'); badperiods = importdata(FileIn); end
            if isstruct(badperiods)
                badperiods = badperiods.data;
            else
                badperiods = [];
            end
        end

        function Tiff2Movie(obj, sourcedir)
            loc = obj.locations;
            if ~exist('sourcedir','var') || isempty(sourcedir)
                sourcedir = fullfiletol(loc.subject_datapath,loc.orig_trials);
            end
            files = obj.filelist;

            parfor i_f = 1:numel(files)
                thisfilename = fullfiletol(sourcedir,files(i_f).name);
                thismatfilename = fullfiletol(sourcedir,[extractBefore(files(i_f).name,'.'),'.mat']);
                if exist(thismatfilename,'file')
                    warning('Output file already found: skipped')
                    continue
                end
                
                thismovie = Movie(thisfilename);

                % make sure the framerate is specified
                if isempty(thismovie.fs)
                    thismovie = thismovie.setFramerateHz(obj.framerate);
                end

                % retrieve badperiods for this trial
                idx = obj.badperiods(:,1) == i_f;
                thismovie.badperiods = obj.badperiods(idx,:);
                
                thismovie.save(sourcedir,'mat')
            end
        end
        
        function localCorr = retrieve_localcorr_maps(obj, input_folder, interval)
            % initialize vars
            if ~exist('input_folder','var') || isempty(input_folder)
                input_folder = obj.locations.orig_trials;
            end
            if ~exist('interval','var') || isempty(interval)
                interval = 1:obj.getNFrames;
            end
            imHeight = obj.singletrial_meta{1}.h;
            imWidth = obj.singletrial_meta{1}.w;
            nTrials = obj.getNTrials;
            localCorr = nan(imHeight,imWidth,nTrials);

            % move to input folder
            cd(fullfiletol(obj.locations.subject_datapath,input_folder))

            % make a list of names for all the files to be loaded
            filenames = {};
            for i_file = 1:nTrials
                filenames{end+1} = find_daughter_file(obj.filelist(i_file).name,'mat');
            end

            parfor i_file = 1:numel(filenames)
                disp(filenames{i_file})
                % assumes that filenames refer to physiological Movie's
                m = robust_io('load',filenames{i_file},'movie').movie;
                M = m.stack(:,:,interval);
                localCorr(:,:,i_file) = computeLocalCorrelationMap(M);
            end
        end

        function saveTrialAvgs(obj, thispath)
            arguments
                obj 
                thispath char = ''
            end

            % init vars
            outpath = 'trialavgs';

            images = obj.retrieve_trial_anatomies(thispath);
        
            % save as a tif movie
            Movie(images).save(outpath,'tif');

            % save individually
            for i_img = 1:size(images,3)
                thisimg = Movie(images(:,:,i_img));
                fname = obj.filelist(i_img).name;
                if endsWith(fname,'.tif'); fname = extractBefore(fname,'.tif'); end
                thisimg.path.fname = fname;
                thisimg.save(fullfiletol(outpath),'tif');
            end
        end

        % Method to retrieve average projections from physiological data.
        % You can manually store the result into obj.anatomy_imgs for 
        % future reference
        function anatomy = retrieve_trial_anatomies(obj,input_folder)
            arguments
                obj 
                input_folder = obj.locations.orig_trials
            end
            ntrials = obj.getNTrials;
            anatomy = nan(height(obj.reference_img),width(obj.reference_img),ntrials);
            
            % move to input folder
            cd(fullfiletol(obj.locations.subject_datapath,input_folder))

            % make a list of names for all the files to be loaded
            filenames = {};
            for i_file = 1:ntrials
                filenames{end+1} = find_daughter_file(obj.filelist(i_file).name,'mat');
            end

            % open each file and extract anatomy image
            parfor i_file = 1:ntrials
                disp(filenames{i_file})
                
                % assumes that filenames refer to physiological Movie's
                m = robust_io('load',filenames{i_file},'movie').movie;
                anatomy(:,:,i_file) = m.timeavg;
            end

        end

        function reference_image = retrieve_ref_img(obj,input_folder, do_save)
            arguments
                obj
                input_folder char
                do_save logical = false
            end

            % specify output location
            output_folder = fullfiletol(obj.locations.subject_datapath,'references',input_folder);

            % move to input folder
            orig_path = pwd;
            cd(fullfiletol(obj.locations.subject_datapath,input_folder))

            reftrialnum = obj.reference_img_meta.Trial_relativenum;
            filename = find_daughter_file(obj.filelist(reftrialnum).name,'mat');

            % extract average of selected frame range
            frame_range = obj.reference_img_meta.Frame_range(1) : ...
                obj.reference_img_meta.Frame_range(2);
            snip = Snippet(filename,frame_range);
            snip.path.fname = getFileNameSpecs(filename).fname; % to have a nice name when saving, if applicable
            reference_image = snip.timeavg;
            
            % return to original folder
            cd(fullfiletol(orig_path))

            if ~do_save; return; end
            
            % Saving as files
            disp(['Saving to ... ', output_folder])
            if ~exist(output_folder,'dir'); mkdir(output_folder); end
            snip.path.fname = ['trial_', num2str(getFileNameSpecs(snip.path.fname).trial_num)];
            snip.save(output_folder,'mat',snip.path.fname)       % Saving the whole Snippet
            FileOut = fullfile(output_folder,[snip.path.fname, '.tif']);    % Saving the image as Tiff
            if exist(FileOut,'file'); delete(FileOut); end            
            saveastiff(reference_image,FileOut)
        end

        function snip = retrieveSnippetfromCoords(obj, filenamejson)
            arguments
                obj
                filenamejson char
        %         Expected json syntax:
        %         {
        %             "Name": "description",
        %             "Trial_relativenum": 35,
        %             "Frame_range": [1050, 1300]
        %         }
            end
            
            meta = readJson(filenamejson);
            
            fname = getFileNameSpecs(obj.filelist(meta.Trial_relativenum).name).fname;
            fname = find_daughter_file(fname,'mat');
            interval = meta.Frame_range(1):meta.Frame_range(2);
            snip = Snippet(fname,interval);
            
            snip.description = meta.Name;
            snip.metadata = meta;
        end

        % Method to log a new state (update the log and currentstate properties)
        function obj = update_currentstate(obj,charv)            
            % append a new voice to log
            try
                obj.log{end+1} = strjoin(charv); % in case of string arrays
            catch
                obj.log{end+1} = charv;
            end
            % make log vertical for easier visualization (only needs to happen once)
            if numel(obj.log)==2; obj.log = transpose(obj.log); end
        
            % change current state
            obj.currentstate = obj.log{end};
        end

        %% Getter methods

        function n = getNTrials(obj)
            n = numel(obj.filelist);
        end

        function L = getNFrames(obj)
            L = obj.singletrial_meta{1}.nfr;
        end

        % Method to set list of files
        function filenames = getFileNames(obj)
            filenames = {obj.filelist(:).name}';
        end
        
        % Method to load the Subject's DATA file (from the Subject's main
        % directory)
        function loadDATAfile(obj)
            loc = obj.locations;

            files = dir(fullfiletol(loc.subject_datapath,'*_DATA*.mat'));
            switch numel(files)
                case 0
                    error('DATA file not found.')
                case 1
                otherwise
                    disp('Multiple DATA files found: ')
                    for i = 1:numel(files); disp(files(i).name); end
                    disp('Last file was loaded.')
            end
            FileIn = fullfiletol(files(1).folder,files(end).name);
            data = robust_io('load',FileIn).data;
            disp('... DATA LOADED.')
        end

        function [fileList, numFiles] = getFileListInSubfolder(obj, folder)
            cd(fullfiletol(obj.locations.subject_datapath))
            
            % retrieve file list
            fileList = {};
            for i = 1:numel(obj.filelist)
                fileList{end+1} = find_daughter_file(...
                    fullfiletol(folder,obj.filelist(i).name),'mat');
            end

            % Determine the number of files
            numFiles = length(fileList);
        end

        %% Functional methods

        % Method to update single-trial ROI definition files with
        % debubbling and center-surround calculations for dF/F calculation.
        % THESE METHODS ACTUALLY NEED VALIDATION
        % function calculate_dFoverF(obj)
        %     original_dir = pwd;
        %     loc = obj.locations;
        %     cd(fullfiletol(loc.subject_datapath,loc.defRois))
        % 
        %     files = dir('*.mat');
        %     for i_f = 1:numel(files)
        %         fname = files(i_f).name;
        %         disp(fname)
        %         load(fname);
        % 
        %         % isolate relevant bad periods
        %         idx = obj.badperiods(:,1) == i_f;
        %         thisbadperiods = obj.badperiods(idx,:);
        % 
        %         % eliminate bad periods from original traces (in case they
        %         % hadn't been already
        %         v = replace_badperiods_with_nans(plane{1}.timetraces,thisbadperiods);
        %         plane{1}.timetraces = v;
        % 
        %         % desurround only
        %         v = desurround2(plane,[],0,1,thisbadperiods);
        %         plane{1}.timetracesdes = v;
        % 
        %         % debubbled desurround
        %         v = debubbled_desurround(plane,num2str(i_f),thisbadperiods);
        %         plane{1}.timetracesdebdes = v;
        % 
        %         save(fname,'plane','-v7.3');
        %     end
        % 
        %     cd(original_dir)
        % end

        function obj = defineTraces(obj)

            tracestmp = ActivityTraces(obj); % construct traces

            % return
            cd(obj.locations.subject_datapath)
            obj.traces = tracestmp;
        end

        % Method to populate DATA structure
        function obj = extract_data(obj)
            disp('-- Extraction of neural activity data --')
            disp(' '); disp(' ')
            
            % inizialize vars
            loc = obj.locations;

            % change directory to the subject's main
            cd(loc.subject_datapath)

            % source single trial data files
            files = dir(fullfiletol(loc.defRois,'*.mat'));

            % identify common units
            disp('Identifying common units.')
            [is_common, ~] = findCommonUnits(fullfiletol(loc.defRois,{files(:).name}),.9,.9);
            disp(' ')
            
            % load sample data for initialization
            disp('Loading sample data for initialization from:')
            FileIn = fullfiletol(loc.defRois,files(1).name);
            disp(FileIn); data = robust_io('load',FileIn).data;
            
            % initialize vars
            % # TODO: for each field, sanity-check if value is the same as in
            % the corresponding Subject attribute
            disp(' ');
            data.L = size(plane{1}.timetraces,1); 
            data.N = numel(is_common);
            [data.traces, data.tracesdes, data.tracesdebdes] = ...
                deal( zeros(data.L*data.N,numel(obj.ROIcheckfiles.to_keep)) );
            data.common_units = is_common;
            data.trials = obj.getFileNames;
            [data.trial_num] = deal(zeros(1,numel(obj.ROIcheckfiles.to_keep)));
            [data.stim_type] = deal({});
            data.meta = plane{1}.meta;
            data.localCorrelations = zeros(size(plane{1}.localCorrelations,1), ...
                size(plane{1}.localCorrelations,2),numel(obj.ROIcheckfiles.to_keep));

            %ROI_map of only common units
            ROI_map = plane{1}.ROI_map;
            isnt_common = setdiff(unique(ROI_map), is_common);
            idx = ismember(ROI_map,isnt_common);
            ROI_map(idx) = 0;
            data.ROI_map_common = ROI_map;
            
            disp(strcat('Period [frames]: ',num2str(data.L)));
            disp(strcat('Period [seconds]: ',num2str(floor(data.L/data.meta.framerate))));
            disp(strcat('Number of common neurons: ',num2str(data.N)));
            disp(strcat('Number of trials: ',num2str(numel(data.trials))));
            
            disp(' ');disp(' ')
            
            for i_f = 1:numel(files)
                filename = fullfiletol(loc.defRois,files(i_f).name);
                disp(strcat('Loading: ',strrep(filename,'\','\\'),' ...'))
                if ~exist(filename,'file'); warning(strcat(strrep(filename,'\','\\'),' not found! - skipped.')); continue; end
                plane = robust_io('load',filename).plane;
                if isempty(plane); error('''plane'' variable not found!'); end

                % original dF/F timetraces
                common_traces = extractCommonTimetracesFromDefROIsFile(plane{1}.timetraces,is_common,data.L);
                data.traces(:,i_f) = common_traces(:);

                % desurrounded dF/F timetraces
                common_traces = extractCommonTimetracesFromDefROIsFile(plane{1}.timetracesdes,is_common,data.L);
                data.tracesdes(:,i_f) = common_traces(:);

                % debubbled-desurrounded dF/F timetraces
                common_traces = extractCommonTimetracesFromDefROIsFile(plane{1}.timetracesdebdes,is_common,data.L);
                data.tracesdebdes(:,i_f) = common_traces(:);

                try
                    data.trial_num(i_f) = getFileNameSpecs(filename).trial_num;
                catch ME
                    data.trial_num(i_f) = str2double(getFileNameSpecs(filename).trial_num);
                end
                data.stim_type{i_f} = obj.stim_series.stimulus{i_f};
                data.localCorrelations(:,:,i_f) = plane{1}.localCorrelations;

                % save fspecs to plane{1}.meta
                plane{1}.meta.fspecs = getFileNameSpecs(filename);
                s.plane = plane;
                robust_io('save',filename,s,'-v7.3');
            end
            
            
            % sort data by trial number
            [data.trial_num, idx] = sort(data.trial_num);
            data.traces = data.traces(:,idx);
            data.stim_type = data.stim_type(idx);
            data.localCorrelations = data.localCorrelations(:,:,idx);
            [data.idx_by_stim_type, ~] = sortbyStimType(data);
            
            % trace quality check (manual 'bad cell' calling)
            % data = checkTracesMan(data);
            
            %  denoise cell traces
            [data.tracesdn, fcutoff] = denoiseCellTraceData(data,0,obj.badperiods,'traces'); % butterworth filt 4 poles
            [data.tracesdesdn, fcutoff] = denoiseCellTraceData(data,0,obj.badperiods,'tracesdes');
            [data.tracesdebdesdn, fcutoff] = denoiseCellTraceData(data,0,obj.badperiods,'tracesdebdes');
            
            [stim_on,stim_off,f0_window,response_window] = obj.set_stimtrig_timestamps('seconds',1);
            
            % downsample denoised traces
            data.meta.downsample = floor(fcutoff*data.meta.framerate/2);
            data.tracesdns = downsample(traceFormat(data.tracesdn,data.L),data.meta.downsample);
            data.Ldns = size(data.tracesdns,1);
            data.tracesdns = traceFormat(data.tracesdns);
            
            data.meta.seriesid = loc.subject_ID;
            data.stim_on_sec = stim_on;
            data.stim_off_sec = stim_off;
            data.f0_window = f0_window;
            data.response_window = response_window;
            
            % select anatomical regions
            data = selectAnatRegions(data,false);
            
            disp('... DONE!'); disp(' ')
            data
            
            % save to Subject file
            obj.data = data;
            obj.save2mat()

            % save again to DATA file
            FileOut = strcat(loc.subject_ID,'_DATA.mat');
            b = prompt_overwrite(FileOut);
            if ~b; return; end
            fprintf(strcat('Saving: ',strrep(FileOut,'\','\\'),' ...'))
            s.data = data;
            robust_io('save',FileOut,s)
            fprintf(' DONE!\n')
        end

        % Method to save object to MAT file
        function save2mat(obj,auto)
            arguments
                obj 
                auto logical = false
            end

            loc = obj.locations;
            fpath = loc.subject_datapath;
            fname = strcat(obj.name,'.mat');
            FileOut = fullfiletol(fpath,fname);
            if auto; b = true; else; b = prompt_overwrite(FileOut); end
            if ~b; return; end
            % eval(strcat(obj.name,'=obj;'));
            s.(obj.name) = obj;
            robust_io('save',FileOut,s,'-v7.3');
        end

        % Method to save a checkfiles file to the subject_datapath location
        function saveROIcheckfiles2mat(obj,auto)            % DELETE this function as soon as it becomes obsolete
            if ~exist('auto','var'); auto = false; end
            loc = obj.locations;
            fname = [loc.subject_ID,'_check.mat'];
            FileOut = fullfiletol(loc.subject_datapath,fname);
            if auto; b = true; else; b = prompt_overwrite(FileOut); end
            if ~b; return; end
            disp(['Saving: ',fname,' ... '])
            s.is_complete = obj.ROIcheckfiles.is_complete;
            s.to_keep = obj.ROIcheckfiles.to_keep;
            s.to_discard = obj.ROIcheckfiles.to_discard;
            robust_io('save',FileOut,s,'-mat')
        end
    
        % Method to reload the Subject MAT file
        function obj = reload(obj)
            loc = obj.locations;
            fname = [obj.name,'.mat'];
            FileIn = fullfiletol(loc.subject_datapath,fname);
            check_exist(FileIn);
            fprintf('Loading: %s',FileIn)
            obj = robust_io('load',FileIn).fish1;
        end
    end
end

%% small accessory functions

function common_traces = extractCommonTimetracesFromDefROIsFile(timetraces,is_common,L)
    common_traces = timetraces(:,is_common(ismember( ...
        is_common,[1:size(timetraces,2)])));
    if size(common_traces,1) > L
        warning(strcat('Period longer than expected: cut to ',num2str(L),' frames.'))
        common_traces = common_traces(1:L,:);
    elseif size(common_traces,1) < L
        warning(strcat('Period of ',num2str(size(common_traces,1)), ...
            ' frames is too short! file skipped.'))
    end
end

