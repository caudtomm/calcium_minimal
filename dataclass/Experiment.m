classdef Experiment
    properties
        name string = ""
        subjectTab table
        subject cell = {}
        locations Locations = Locations
        summaryTab table = table()
    end
    
    methods 
        function obj = Experiment(tabpath,sheet)
            arguments
                tabpath char {isfile}
                sheet char = ''
            end
            obj.subjectTab = loadContentTab(tabpath,sheet);
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

        function obj = processSubjects(obj, subjectlist)
            
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