classdef Subject
    properties
        locations
        id
        name
        log
        imaging_date {mustBeNumeric}
        group
        data
        odor_delay {mustBeNonnegative}
        notes
        filelist
        reference_img {mustBeNumeric}
        reference_img_meta
        anatomy_imgs
        trialnum
        min_trialnum
        stim_series
        ch2stimtype_map
        framerate
        scanimage_metadata
        badperiods
        currentstate
        ROIcheckfiles
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
            [obj.reference_img, obj.reference_img_meta] = obj.loadReferenceImg(obj.locations.references.raw);
            obj.anatomy_imgs = [];
            obj.badperiods = obj.loadBadperiods();
            obj.log = {};
            obj = update_currentstate(obj, 'newly constructed');
            obj.ROIcheckfiles = struct('to_keep',[],'to_discard',[],'is_complete',[]);

            obj.save2mat();
        end
        
        %% Setter methods

        % Method to set list of files
        function [files,trialnum,min_trialnum] = setFileList(obj)
            loc = obj.locations;
            relative_pathin = fullfile(loc.subject_ID, loc.rawtrials);

            fprintf(strcat('\nSource folder: ', strrep(relative_pathin, '\','\\'),'.\n\n'))

            files = dir(fullfile(loc.general_datapath, relative_pathin, ...
                strcat('*.',loc.datafile_ext)));

            % order by trial number
            n = [];
            for i_f = 1:numel(files)
                file = files(i_f);
                fspecs = getFileNameSpecs(file.name);
                n = [n;fspecs.trial_num];
            end
            min_trialnum = min(n)-1; n=n-min_trialnum;
            [~,idx] = sort(n); trialnum = n(idx);
            files = files(idx);
        end
        
        % Method to read series of stimuli (output file of LabVIEW odor 
        % system control script)
        function [stim_series, ch2stimtype_map] = setStimSeries(obj)
            loc = obj.locations;
            fname = strcat(loc.subject_ID,'_',num2str(obj.min_trialnum+1, '%05.f'));
            FileIn = fullfile(loc.subject_datapath,loc.rawtrials,fname);
            [stim_series, ch2stimtype, ch2stimtype_map] = read_stim_series(FileIn);
            stim_series = horzcat(array2table(stim_series, ...
                'VariableNames',{'trialnum' 'odor_channel' 'frame_onset','frame_offset'}), ...
                cell2table(ch2stimtype,'VariableNames',{'stimulus'}));
        end

        % Method to read imaging metadata from scanimage-a
        function [framerate, scanimage_metadata] = setScanimageMetadata(obj) % FROM THE FIRST FILE
            FileIn = fullfile(obj.filelist(1).folder,obj.filelist(1).name);
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

        % Method to load and set the calcium imaging reference image (1 for
        % each plane, concatenated in Z)
        function [img,img_meta] = loadReferenceImg(obj,subdirectory)
            disp(['Loading reference image from: ', subdirectory]);

            [img,img_meta] = deal([]);
            loc = obj.locations;
            PathIn = fullfile(loc.subject_datapath, subdirectory);
            
            % load reference image
            fname = strcat('*.',loc.datafile_ext);
            files = dir(fullfile(PathIn,fname));
            switch numel(files)
                case 0
                    warn_notfound('raw reference image'); return
                case 1
                otherwise
                    warning('more than one reference image available. selected: %s', ...
                        files(1).name)
            end
            img = double(loadTiffStack(char(fullfile(PathIn,files(1).name))));
            
            % load metadata
            fname = 'metadata.json';
            FileIn = fullfile(PathIn,fname);
            if ~exist(FileIn,'file'); warn_notfound(fname); return; end
            jsonStr = fileread(FileIn);
            try
                img_meta = jsondecode(jsonStr);
                showcontent(img_meta)
            catch
                warning('couldn''t read: %s', fname)
            end
        end
        
        % Method to load BADPERIOD coordinates
        function badperiods = loadBadperiods(obj)
            loc = obj.locations;
            FileIn = fullfile(loc.subject_datapath,loc.rawtrials,[loc.subject_ID,'_badperiods.csv']);
            badperiods.data = [];
            if exist(FileIn,'file'); badperiods = importdata(FileIn); end
            if isstruct(badperiods)
                badperiods = badperiods.data;
            else
                badperiods = [];
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
            FileIn = fullfile(obj.locations.subject_datapath,obj.locations.defRois,[fname,'_defROIs.mat']);
            load(FileIn);
            p_ann = plane{1}.ROI_map;
        end
        
        % Method to retrieve average projections from physiological data.
        % You can manually store the result into obj.anatomy_imgs for 
        % future reference
        function anatomy = retrieve_trial_anatomies(obj,input_folder)
            loc = obj.locations;
            ntrials = numel(obj.trialnum);

            % move to input folder
            cd(fullfile(obj.locations.subject_datapath,input_folder))

            % make a list of names for all the files to be loaded
            filenames = {};
            for i_file = 1:numel(obj.filelist)
                filenames{end+1} = find_daughter_file(obj.filelist(i_file).name,loc.datafile_ext);
            end
            
            % open each file and extract anatomy image
            anatomy = nan(height(obj.reference_img),width(obj.reference_img),ntrials);
            for i_file = 1:ntrials
                % assumes that filenames refer to physiological stacks
                movie = double(loadTiffStack(filenames{i_file}));
                anatomy(:,:,i_file) = nanmean(movie,3);
            end

            % plot anatomies vs reference image
%             obj.visualize_anatomy_physFOV;
        end

        % Method to log a new state (update the log and currentstate properties)
        function obj = update_currentstate(obj,charv)
            % append a new voice to log
            obj.log{end+1} = charv;
            % make log vertical for easier visualization (only needs to happen once)
            if numel(obj.log)==2; obj.log = transpose(obj.log); end
        
            % change current state
            obj.currentstate = obj.log{end};
        end

        % Method to update the object property 'ROIcheckfiles'
        function obj = update_ROIcheckfiles(obj,filename,auto)
            thisiscomplete = true; if ~auto; thisiscomplete = askBoolean('Annotation complete? [y,n] '); end
            if thisiscomplete && ~any(strcmp(obj.ROIcheckfiles.is_complete,filename))
                obj.ROIcheckfiles.is_complete{end+1} = filename;
            end

            isgood = true; if ~auto; isgood = askBoolean('Trial good? [y,n] '); end
            if isgood && ~any(strcmp(obj.ROIcheckfiles.to_keep,filename))
                obj.ROIcheckfiles.to_keep{end+1} = filename;
            elseif ~isgood && ~any(strcmp(obj.ROIcheckfiles.to_discard,filename))
                obj.ROIcheckfiles.to_discard{end+1} = filename;
            end
        end


        %% Getter methods

        % Method to set list of files
        function filenames = getFileNames(obj)
            filenames = {obj.filelist(:).name}';
        end
        
        % Method to load the Subject's DATA file (from the Subject's main
        % directory)
        function loadDATAfile(obj)
            loc = obj.locations;

            files = dir(fullfile(loc.subject_datapath,'*_DATA*.mat'));
            switch numel(files)
                case 0
                    error('DATA file not found.')
                case 1
                otherwise
                    disp('Multiple DATA files found: ')
                    for i = 1:numel(files); disp(files(i).name); end
                    disp('Last file was loaded.')
            end
            FileIn = fullfile(files(1).folder,files(end).name);
            load(FileIn)
            disp('... DATA LOADED.')
        end

        %% Functional methods

        % Method to register each imaging plane to a common
        % reference across files
        function [obj,reg_stack,reg_histeqstack,allstats] = registration(obj, method)
            loc = obj.locations;
            loc_raw = fullfile(loc.subject_datapath,loc.rawtrials);
            loc_rigidreg_raw = fullfile(loc.subject_datapath,loc.rawtrials_rigidreg);
            loc_histeq = fullfile(loc.subject_datapath,loc.histeqtrials);

            if isempty(obj.reference_img); error('No reference image found!'); end

            allstats = {};

            for i_file = 1:numel(obj.filelist)
                FilenameRAW = fullfile(loc_raw,obj.filelist(i_file).name);
                FilenameHISTEQ = fullfile(loc_histeq,obj.filelist(i_file).name);
                fspecs = getFileNameSpecs(char(FilenameRAW));

                thisbadperiods = [];
                if ~isempty(obj.badperiods); thisbadperiods = obj.badperiods(obj.badperiods(:,1)==fspecs.trial_num,:); end
                
                switch method
                    case 'rigid'
                        [reg_stack,reg_histeqstack,stats] = ...
                            registration_rigid_singlestack( ...
                            FilenameRAW, ...
                            FilenameHISTEQ, ...
                            obj.reference_img, ...
                            thisbadperiods, ...
                            obj.framerate);
                        allstats{end+1} = stats;
                    case 'opticflow'
                        cd(loc_rigidreg_raw)
                        FilenameRAW_rigidreg = find_daughter_file(obj.filelist(i_file).name,loc.datafile_ext);
                        registration_opticflow_singlestack( ...
                            FilenameRAW_rigidreg, ...
                            obj.reference_img);
                    case 'bunwarpj' % CREATE FUNCTION FROM STACK_DEWOBBLER
                        cd(loc_rigidreg_raw)
                        FilenameRAW_rigidreg = find_daughter_file(obj.filelist(i_file).name,loc.datafile_ext);
                        registration_bunwarpj_singlestack();
                    otherwise
                        error('registration method unknown')
                end
            end

            switch method
                case 'opticflow'
                    cd(loc_rigidreg_raw)
                    registration_opticflow_moveresults(loc,true)
                case 'bunwarpj'
                    % script to concatenate and move dewobbled stacks
                    
            end

            obj = obj.update_currentstate(['registration complete - method: ',method]);
        end 
        
        % Method to visualize anatomical projections from all trials
        function visualize_anatomy_physFOV(obj)
            % based on the output of obj.retrieve_trial_anatomies()

            anatomy = obj.anatomy_imgs;
            ntrials = numel(obj.trialnum);

            % figure 1
            figure;
            for i_trial = 1:ntrials
                subplot(6,6,i_trial); imagesc(anatomy(:,:,i_trial));
                title(num2str(i_trial))
            end
            subplot(6,6,36); imagesc(obj.reference_img); title('reference')
        end

        % Method containing the FULL PROTOCOL to define ROIs
        function [obj] = define_ROIs(obj,input_folder,mode)
            % move to input folder
            cd(fullfile(obj.locations.subject_datapath,input_folder))

            disp('-- ROI selection --');disp('')
            
            if ~exist('mode','var') || isempty(mode)
                mode = prompt_chooseString('Specify selection mode', {'manual','automatic'});
            end
            obj.define_ROIs_oneway(mode);

            
        end

        % Method to define ROIs manually through a list of files
        function obj = define_ROIs_oneway(obj,mode)
            loc = obj.locations;
            fs = obj.framerate;
            
            % implement mode-specific options
            switch mode
                case 'automatic'
                    auto = true;
                case 'manual'
                    auto = false;
                otherwise
                    error('ROI selection mode unknown.')
            end

            % initialize looping vars
            partial_annotation = obj.load_partial_annotation();
            go_on = true;

            for i_file = 1:numel(obj.filelist)
                %skip if already complete
                if ~auto && any(ismember(obj.filelist(i_file).name,obj.ROIcheckfiles.is_complete)); continue; end
                
                filename = find_daughter_file(obj.filelist(i_file).name,loc.datafile_ext);

                [~,~,f0_window,response_window] = obj.set_stimtrig_timestamps('frames',i_file);
                
                plane = defROIs2(filename, ...
                    obj.scanimage_metadata, ...
                    partial_annotation, ...
                    auto, ... % do_prune
                    f0_window, ...
                    response_window, ...
                    auto); % auto
                
                % save
                fspecs = getFileNameSpecs(obj.filelist(i_file).name);
                FileOut_path = fullfile(loc.subject_datapath,loc.defRois);
                mkdir(FileOut_path)
                FileOut = fullfile(FileOut_path,[fspecs.fname,'_defROIs.mat']);
                if auto; b = true; else; b = prompt_overwrite(FileOut); end
                if b; save(FileOut,'plane','-v7.3'); end

                % update vars and attributes
                partial_annotation = plane{1}.ROI_map;
                obj = obj.update_ROIcheckfiles(obj.filelist(i_file).name,auto);
                
                % save in itinere (lose the account of maximum one trial!)
                obj.save2mat(true);
                
                % continue to next trial?
                if ~auto; go_on = askBoolean('Wanna continue to the next trial? [y,n] '); end
                if ~go_on; break; end
   
            end
        end

        % Method to update single-trial ROI definition files with
        % debubbling and center-surround calculations for dF/F calculation.
        % THESE METHODS ACTUALLY NEED VALIDATION
        function calculate_dFoverF(obj)
            original_dir = pwd;
            loc = obj.locations;
            cd(fullfile(loc.subject_datapath,loc.defRois))

            files = dir('*.mat');
            for i_f = 1:numel(files)
                fname = files(i_f).name;
                disp(fname)
                load(fname);

                % isolate relevant bad periods
                idx = obj.badperiods(:,1) == i_f;
                thisbadperiods = obj.badperiods(idx,:);
                
                % eliminate bad periods from original traces (in case they
                % hadn't been already
                v = replace_badperiods_with_nans(plane{1}.timetraces,thisbadperiods);
                plane{1}.timetraces = v;

                % desurround only
                v = desurround2(plane,[],0,1,thisbadperiods);
                plane{1}.timetracesdes = v;

                % debubbled desurround
                v = debubbled_desurround(plane,num2str(i_f),thisbadperiods);
                plane{1}.timetracesdebdes = v;
            
                save(fname,'plane','-v7.3');
            end

            cd(original_dir)
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
            files = dir(fullfile(loc.defRois,'*.mat'));

            % identify common units
            disp('Identifying common units.')
            [is_common, ~] = findCommonUnits(fullfile(loc.defRois,{files(:).name}),.9,.9);
            disp(' ')
            
            % load sample data for initialization
            disp('Loading sample data for initialization from:')
            FileIn = fullfile(loc.defRois,files(1).name);
            disp(FileIn); load(FileIn)
            
            % initialize vars
            % TODO: for each field, sanity-check if value is the same as in
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
                filename = fullfile(loc.defRois,files(i_f).name);
                disp(strcat('Loading: ',strrep(filename,'\','\\'),' ...'))
                if ~exist(filename,'file'); warning(strcat(strrep(filename,'\','\\'),' not found! - skipped.')); continue; end
                load(filename)
                if ~exist('plane','var'); error('''plane'' variable not found!'); end
                
%                 common_traces = plane{1}.timetraces(:,is_common(ismember( ...
%                     is_common,[1:size(plane{1}.timetraces,2)])));
%                 if size(common_traces,1) > data.L
%                     warning(strcat('Period longer than expected: cut to ',num2str(data.L),' frames.'))
%                     common_traces = common_traces(1:data.L,:);
%                 elseif size(common_traces,1) < data.L
%                     warning(strcat('Period of ',num2str(size(common_traces,1)), ...
%                         ' frames is too short! file skipped.'))
%                 end

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
                save(filename,'plane','-v7.3');
            end
            
            
            % sort data by trial number
            [data.trial_num, idx] = sort(data.trial_num);
            data.traces = data.traces(:,idx);
            data.stim_type = data.stim_type(idx);
            data.localCorrelations = data.localCorrelations(:,:,idx);
            [data.idx_by_stim_type, ~] = sortbyStimType(data);
            
            % trace quality check (manual 'bad cell' calling)
%             data = checkTracesMan(data);
            
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
            save(FileOut,'data')
            fprintf(' DONE!\n')
        end

        % Method to save object to MAT file
        function save2mat(obj,auto)
            if ~exist('auto','var'); auto = false; end
            loc = obj.locations;
            fpath = loc.subject_datapath;
            fname = strcat(obj.name,'.mat');
            FileOut = fullfile(fpath,fname);
            if auto; b = true; else; b = prompt_overwrite(FileOut); end
            if ~b; return; end
            eval(strcat(obj.name,'=obj;'));
            eval(strcat('save(FileOut,''',obj.name,''')'));
        end
        
        % Method to save a checkfiles file to the subject_datapath location
        function saveROIcheckfiles2mat(obj,auto)            % DELETE this function as soon as it becomes obsolete
            if ~exist('auto','var'); auto = false; end
            loc = obj.locations;
            fname = [loc.subject_ID,'_check.mat'];
            FileOut = fullfile(loc.subject_datapath,fname);
            if auto; b = true; else; b = prompt_overwrite(FileOut); end
            if ~b; return; end
            disp(['Saving: ',fname,' ... '])
            is_complete = obj.ROIcheckfiles.is_complete;
            to_keep = obj.ROIcheckfiles.to_keep;
            to_discard = obj.ROIcheckfiles.to_discard;
            save(FileOut,'is_complete',"to_keep","to_discard",'-mat')
        end
    
        % Method to reload the Subject MAT file
        function obj = reload(obj)
            loc = obj.locations;
            fname = [obj.name,'.mat'];
            FileIn = fullfile(loc.subject_datapath,fname);
            check_exist(FileIn);
            fprintf('Loading: %s',FileIn)
            load(FileIn)
        end
    end
end

%% small accessory functions

function warn_notfound(filename)
    warning('File not found: %s', filename);
end

function showcontent(obj)
    obj
end

function b = prompt_overwrite(filename)
    b = true;
    if ~exist(filename,'file'); return; end
    str = ['Pre-existing file found: \n',char(filename),'\nOverwrite? [y/n]'];
    b = askBoolean(str);
    if ~b; warning('Aborted overwrite.'); end
end

function daughterfile = find_daughter_file(basename,obligate_ext)
    % search is carried out in the current directory
    if ~isempty(obligate_ext)
        files = dir(['*.',obligate_ext]);
    else 
        files = dir();
    end
    files = {files(:).name}';

    % remove extension (if applicable)
    if contains(basename,'.'); basename = extractBefore(basename,'.'); end

    % search daughters
    idx = contains(files, basename);
    daughterfile = files{find(idx,1)};
    if sum(idx)>1; warning('%s daugther files found: used %s',num2str(sum(idx)),daughterfile); end
end

function check_exist(FileIn)
    if ~exist("FileIn")
        str = ['Not found: ', FileIn];
        error(str)
    end
end

function showtrial(src,event,i_trial)
    switch event.Key
        case {'downarrow','rightarrow'}
            disp(i_trial)
        case {'uparrow','leftarrow'}
            disp('bye')
    end
end

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

