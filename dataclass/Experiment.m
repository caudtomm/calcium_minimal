classdef Experiment
    properties
        name string = ""
        tracesfolder char
        record
        subjectTab table
        subject cell = {}
        currentsubject Subject
        traces cell = {}
        locations Locations = Locations
        summaryTab table = table()
    end
    
    methods (Access = private)
        function p = preprocessingHead(obj)
            % startup preprocessing
            p = Preprocessing;
            p.autosave = true;
            p.pl = false;
    
            % load subject into Preprocessing
            p.sj = obj.currentsubject;
    
            % % make sure the subject has the correct reference image
            % p = p.updateSubject(p.sj.locations.rawtrials);
        end
    end

    methods 
        function obj = Experiment(tabpath,tracesfolder,sheet)
            arguments
                tabpath char {isfile}
                tracesfolder char = 'traces\top1biasedIC'
                sheet char = ''
            end
            obj.record.tabpath = tabpath;
            obj.record.sheet = sheet;
            obj.subjectTab = loadContentTab(tabpath,sheet);
            obj.name = string(sheet);
            obj.tracesfolder = tracesfolder;
        end

        function obj = loadSubjects(obj,subjectnames,create,do_overwrite)
            arguments
                obj
                subjectnames cell = obj.subjectTab.name
                create logical = true;
                do_overwrite logical = false; % if true, previously loaded subject are discarded
            end

            % initialize
            loc = obj.locations;
            nsubjects = numel(subjectnames);
            fname = 'fish1.mat';
            if do_overwrite; obj.subject = {}; end

            % execute
            for i_subject = 1:nsubjects
                disp('')

                thissubjectname = subjectnames{i_subject};
                disp(thissubjectname)
                idx = strcmp(obj.subjectTab.name,thissubjectname);

                % create temporary Locations object for the Subject
                thisloc = loc; % copy all the information of the experiment's Locations
                thisloc = thisloc.setSubjectID(thissubjectname);

                % look for pre-existing subject file
                FileIn = fullfiletol(thisloc.subject_datapath,fname);
                if isfile(FileIn)
                    fprintf('found : %s ...',fname)
                    thissubject = robust_io('load',FileIn,'fish1').fish1;
                    fprintf(' loaded.')
                    disp('')
                elseif create
                    fprintf('not found : %s ...',fname)
                    thisgroup = obj.subjectTab.group(idx);
                    thisgroup = thisgroup(1); % for robustness
                    thissubject = Subject('fish1',thisloc,thisgroup);
                    fprintf(' created.')
                    disp('')
                else
                    fprintf('not found : %s ...',fname)
                    fprintf(' skipped.')
                    disp('')
                    continue
                end

                % take local tungsten drive name and data folder into account
                thissubject.locations = ...
                    thissubject.locations.setDrive(obj.locations.drive);
                thissubject.locations = ...
                    thissubject.locations.setDataFolder(obj.locations.datafolder);

                % save to properties                
                subject_not_yet_in_this_experiment = ...
                    isempty(obj.summaryTab) || ...
                    sum(strcmp(obj.summaryTab.name,thissubjectname))==0;

                if subject_not_yet_in_this_experiment
                    obj.subject{end+1} = thissubject;
                    obj.summaryTab(end+1,:) = obj.subjectTab(idx,:);
                else
                    disp('pre-existing: OVERWRITE')

                    idx_preexisting = strcmp({obj.subject.id},thissubjectname);
                    obj.subject{idx_preexisting} = thissubject;

                    idx_preexisting = strcmp(obj.summaryTab.name,thissubjectname);
                    obj.summaryTab(idx_preexisting,:) = obj.subjectTab(idx,:);
                end
                    
            end

            % update current subject
            obj.currentsubject = obj.subject{end};
        end

        function obj = loadSubjectTraces(obj,subjectnames,do_overwrite)
             arguments
                obj
                subjectnames cell = obj.subjectTab.name
                do_overwrite logical = false; % if true, previously loaded traces are discarded
             end
            
            % initialize
            loc = obj.locations;
            nsubjects = numel(subjectnames);
            if do_overwrite; obj.traces = {}; end

            % execute
            for i_subject = 1:nsubjects
                disp('')
                disp('')

                thissubjectname = subjectnames{i_subject};
                disp(thissubjectname)

                % create Locations object for the Subject
                thisloc = loc; % copy all the information of the experiment's Locations
                thisloc = thisloc.setSubjectID(thissubjectname);

                % look for pre-existing traces file
                % fname = [thissubjectname,'_traceslight.mat'];
                % FileIn = fullfiletol(thisloc.subject_datapath,fname);
                files = dir(fullfiletol(thisloc.subject_datapath,obj.tracesfolder,'*traceslight*'));
                idx = arrayfun(@(x) endsWith(x.name,'.mat'),files);
                files = files(idx);
                if ~isempty(files)
                    thisfile = files(1);
                    fname = thisfile.name;
                    FileIn = fullfiletol(thisfile.folder,thisfile.name);
                    fprintf('found : %s ...',fname)
                    thistraces = robust_io('load',FileIn).traces;
                    fprintf(' loaded.')
                    disp(''); disp('')
                else
                    disp('no traces found! skipped!')
                    disp(''); disp('')
                    continue
                end

                % take local tungsten drive name into account
                thistraces.subject_locations = ...
                    thistraces.subject_locations.setDrive(obj.locations.drive);


                % save to properties
                
                % check if already loaded
                idx = numel(obj.traces)+1;
                idx_matching_subjectname = find(cellfun(@(x) ...
                    strcmp(x.subject_locations.subject_ID, thissubjectname), obj.traces),1);
                if ~isempty(idx_matching_subjectname)
                    disp(['found preloaded traces for : ',thissubjectname])
                    disp('OVERWRITING')
                    idx = idx_matching_subjectname;
                end
                % actually save
                obj.traces{idx} = thistraces;
                    
            end

            disp('')
            disp('')
        end
        
        function obj = processSubjects(obj,process_name,docheckifdone)
            arguments
                obj 
                process_name char
                docheckifdone logical = true
            end
            % define complete list of subjects to process
            subjectlist = obj.subjectTab.name;

            for i_sub = 1:numel(subjectlist)

                % check if already done: if so, skip
                if docheckifdone & strcmp(obj.subjectTab.notes{i_sub},'done!')
                    continue
                end

                % load the current subject and move to its data directory
                obj = obj.loadSubjects(subjectlist(i_sub),false,true);
                cd(obj.currentsubject.locations.subject_datapath)
                
                % --------------------------------------
                % VARIABLE SECTION
                %
                p = obj.preprocessingHead();
                validMethods = methods(p);
                
                if ismember(process_name, validMethods)
                    idx_tosubtract = str2double(obj.subjectTab.manPickedICs_conservative_{i_sub});
                    if idx_tosubtract==0; continue; end

                    p = p.(process_name)(idx_tosubtract); % add any input arguments on the right
                else
                    error('process unknown.')
                end
                %
                % --------------------------------------

                obj.subjectTab.notes{i_sub} = 'done!';
                writetable(obj.subjectTab,obj.record.tabpath,'Sheet',obj.record.sheet)
            end
            
        end

        function obj = saveTrialAvgsToFiles(obj)
            % define complete list of subjects to process
            subjectlist = obj.subjectTab.name;

            for i_sub = 1:numel(subjectlist)
                disp(subjectlist(i_sub))
                
                % load the current subject and move to its data directory
                obj = obj.loadSubjects(subjectlist(i_sub),false,true);
                cd(obj.currentsubject.locations.subject_datapath)

                folders = dir(pwd); % includes files, too.
                for i = 3:numel(folders)
                    if ~folders(i).isdir; continue; end

                    obj.currentsubject.saveTrialAvgs(folders(i).name);
                end
            end
        end

        function experiment = convert2BackwardCompatibleStruct(obj, name)
            arguments
                obj 
                name char = 'new_experiment'
            end

            experiment.name = name;
            experiment.summaryTable = obj.subjectTab;

            for i = 1:numel(obj.traces)
                 experiment.series{i} = obj.traces{i}.convert2Series;
            end
        end

        function obj = saveAnatomicalReferencesToFiles(obj)
            % define complete list of subjects to process
            subjectlist = obj.subjectTab.name;

            for i_sub = 1:numel(subjectlist)
                disp(subjectlist(i_sub))

                % load the current subject and move to its data directory
                obj = obj.loadSubjects(subjectlist(i_sub),false,true);
                cd(obj.currentsubject.locations.subject_datapath)

                folders = dir(pwd); % includes files, too.
                for i = 3:numel(folders)
                    if ~folders(i).isdir; continue; end

                    try % some folders don't actually contain trials
                        obj.currentsubject.retrieve_ref_img(folders(i).name,true);
                    catch
                        % just skip
                    end
                end             
                
            end
            
        end

        
    end
end

function contentTab = loadContentTab(tabpath,sheet)
    arguments
        tabpath char {isfile}
        sheet char = ''
    end
    
    if isempty(sheet)
        contentTab = readtable(tabpath);
    else
        contentTab = readtable(tabpath,'Sheet',sheet);
    end
end