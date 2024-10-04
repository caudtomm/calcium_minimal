classdef Experiment
    properties
        name string = ""
        record
        subjectTab table
        subject cell = {}
        currentsubject Subject
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
    
            % make sure the subject has the correct reference image
            p = p.updateSubject(p.sj.locations.rawtrials);
        end
    end

    methods 
        function obj = Experiment(tabpath,sheet)
            arguments
                tabpath char {isfile}
                sheet char = ''
            end
            obj.record.tabpath = tabpath;
            obj.record.sheet = sheet;
            obj.subjectTab = loadContentTab(tabpath,sheet);
            obj.name = string(sheet);
        end

        function obj = loadSubjects(obj,subjectnames,create)
            arguments
                obj
                subjectnames cell = obj.subjectTab.name
                create logical = true;
            end

            % initialize
            loc = obj.locations;
            nsubjects = numel(subjectnames);
            fname = 'fish1.mat';

            % execute
            for i_subject = 1:nsubjects
                disp('')

                thissubjectname = subjectnames{i_subject};
                disp(thissubjectname)
                idx = strcmp(obj.subjectTab.name,thissubjectname);

                % create Locations object for the Subject
                thisloc = loc; % copy all the information of the experiment's Locations
                thisloc = thisloc.setSubjectID(thissubjectname);

                % look for pre-existing subject file
                FileIn = fullfile(thisloc.subject_datapath,fname);
                if isfile(FileIn)
                    fprintf('found : %s ...',fname)
                    thissubject = load(FileIn,'fish1').fish1;
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

                % load the current subject
                obj = obj.loadSubjects(subjectlist(i_sub));
                obj.currentsubject = obj.subject{1};
                
                obj.currentsubject.locations = ...
                    obj.currentsubject.locations.setDrive(obj.locations.drive);
                cdtol(newlocation.subject_datapath)

                % --------------------------------------
                % VARIABLE SECTION
                %
                p = obj.preprocessingHead();
                try
                    p = p.(process_name); % ###
                catch
                    error('process unknown.')
                end
                %
                % --------------------------------------

                obj.subjectTab.notes{i_sub} = 'done!';
                writetable(obj.subjectTab,obj.record.tabpath,'Sheet',obj.record.sheet)

                obj.subject = {};
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