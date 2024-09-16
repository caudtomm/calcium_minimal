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
        
        rawtrials char = 'trials'
        references struct = struct()
        histeqtrials = 'trials_clahe';
        rawtrials_rigidreg = 'reg_stacks_raw';
        histeqtrials_rigidreg = 'reg_stacks_clahe';
        rawtrials_opticflowwarp = 'OFreg_stacks_raw';
        histeqtrials_opticflowwarp = 'OFreg_stacks_clahe';

%         locations.warp.rawtrials_bunwarpj = 'trials_warp_bunwarpj';
        
    end

    methods
        function obj = Locations()
            obj = obj.setGeneralDataPath(...
                fullfile(obj.drive, ...
                obj.datafolder));

            obj.references.raw = fullfile('references','rawtrials');
            obj.references.histeq = fullfile('references','histeqtrials');
        end

        function obj = setDrive(obj, id)
            arguments
                obj
                id char = ''
            end
            
            if ~isempty(id)
                obj.drive = id;
            end % else, keep existing

            obj = obj.setGeneralDataPath(fullfile(obj.drive,obj.datafolder));

        end

        function obj = setGeneralDataPath(obj, id)
            arguments
                obj
                id char = ''
            end

            if ~isempty(id)
                obj.general_datapath = id;
            end % else, keep existing

            obj = obj.setSubjectDataPath(fullfile(obj.general_datapath,obj.subject_ID));

        end
        
        function obj = setSubjectDataPath(obj, id)
            arguments
                obj
                id char
            end

            obj.subject_datapath = id;

        end

        function obj = setSubjectID(obj, id)
            arguments
                obj
                id char = ''
            end

            if ~isempty(id)
                obj.subject_ID = id;
            end % else, keep existing

            obj = obj.setSubjectDataPath(fullfile(obj.general_datapath,obj.subject_ID));

        end

        function obj = setDataFolder(obj, id)
            arguments
                obj
                id char = ''
            end

            if ~isempty(id)
                obj.datafolder = id;
            end % else, keep existing

            obj = obj.setGeneralDataPath(fullfile(obj.drive,obj.datafolder));

        end
    end
end