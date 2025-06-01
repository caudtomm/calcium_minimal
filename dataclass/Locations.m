classdef Locations
    properties (SetAccess = private)
        drive = 'W:'
        general_datapath char
        subject_datapath char
        subject_ID char = ''
        datafolder char = 'scratch\gfriedri\caudtomm\ev_data'
    end

    properties
        datafile_ext char = 'tif'

        traces_src char = ''    % this is the folder with the movies used 
                                % for the extraction of calcium traces:
                                % set it manually, preferrably by pointing
                                % to another specified folder name
                                % property of the same Location object
        
        orig_trials char = 'orig_trials'
        rawtrials char = 'trials'
        references struct = struct()
        histeqtrials = 'trials_clahe';
        rawtrials_rigidreg = 'reg_stacks_raw';
        histeqtrials_rigidreg = 'reg_stacks_clahe';
        rawtrials_rigidreg_fromhisteq = 'reg_stacks_clahe2raw';
        rawtrials_opticflowwarp = 'OFreg_stacks_raw';
        histeqtrials_opticflowwarp = 'OFreg_stacks_clahe';
        rawtrials_opticflowwarp_fromhisteq = 'OFreg_stacks_clahe2raw';
        
        rawtrials_opticflowwarp_fromhisteq_corrected = 'OFreg_stacks_clahe2raw_corr';
        rawtrials_opticflowwarp_fromhisteq_corrected_noIC = 'OFreg_stacks_clahe2raw_IC';
        
    end

    methods (Static)
        function pathOut = tungstenDrive(pathIn)
            arguments
                pathIn char = ''
            end

            persistent pathPersistent;
            pathPersistent = char(pathPersistent);

            if ~isempty(pathIn)
                pathPersistent = pathIn;
            end

            pathOut = pathPersistent;

        end
    end

    methods
        function obj = Locations()
            obj = obj.setGeneralDataPath(...
                fullfiletol(obj.drive, ...
                obj.datafolder));

            obj.references.raw = fullfiletol('references','rawtrials');
            obj.references.histeq = fullfiletol('references','histeqtrials');
        end

        function obj = setDrive(obj, id)
            arguments
                obj
                id char = ''
            end
            
            % This will retain the existing drive, in case the id was
            % empty.
            obj.drive = obj.tungstenDrive(fullfiletol(id));

            obj = obj.setGeneralDataPath(fullfiletol(obj.drive,obj.datafolder));

        end

        function obj = setGeneralDataPath(obj, id)
            arguments
                obj
                id char = ''
            end

            if ~isempty(id)
                obj.general_datapath = fullfiletol(id);
            end % else, keep existing

            obj = obj.setSubjectDataPath(fullfiletol(obj.general_datapath,obj.subject_ID));

        end
        
        function obj = setSubjectDataPath(obj, id)
            arguments
                obj
                id char
            end

            obj.subject_datapath = fullfiletol(id);

        end

        function obj = setSubjectID(obj, id)
            arguments
                obj
                id char = ''
            end

            if ~isempty(id)
                obj.subject_ID = fullfiletol(id);
            end % else, keep existing

            obj = obj.setSubjectDataPath(fullfiletol(obj.general_datapath,obj.subject_ID));

        end

        function obj = setDataFolder(obj, id)
            arguments
                obj
                id char = ''
            end

            if ~isempty(id)
                obj.datafolder = fullfiletol(id);
            end % else, keep existing

            obj = obj.setGeneralDataPath(fullfiletol(obj.drive,obj.datafolder));

        end
    end
end