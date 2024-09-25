classdef ProcessInit
    properties
        description string = ""
        detailed_description string = ""
        processclass char = ''
        original_path char = pwd
        original_fname char = ''
        hash char = ''
        superHash char = ''
        subHash cell = {}
        batch_run logical = false
        operation_precomputed logical = false
        orig_computeHash char = ''
        keep_heavy_res_in_memory logical = true;
        date_firstrun datetime = datetime("today")
        method char = ''
        reference_meta struct
        internal_settings
    end

    methods (Static)
        function newobj = fromStruct(struct_in)
            newobj = ProcessInit();
            for fn = fieldnames(struct_in)'    %enumerat fields
              try
                  newobj.(fn{1}) = struct_in.(fn{1});   %and copy
              catch
                  warning('Could not copy field %s', fn{1});
              end
           end
        end
    end
    methods
        function obj = ProcessInit()
            
        end

        function save(obj,newpath,type) % type : {'mat', json'}
            arguments
                obj ProcessInit
                newpath char = pwd
                type char = 'json'
            end

            fname = ['',obj.hash,'_init'];
            outpath = newpath;

            init = obj;
            switch type
                case 'mat'
                    % save init as mat file
                    FileOut = fullfile(outpath,[fname,'.mat']);
                    b = prompt_overwrite(FileOut);
                    if b
                        save(FileOut,'init');
                    end
                case 'json'
                    % save init as json file
                    FileOut = fullfile(outpath,[fname,'.json']);
                    b = prompt_overwrite(FileOut);
                    if b
                        fid = fopen(FileOut,'w');
                        if fid == -1
                            error('File could not be opened, check permissions.');
                        end
                        formattedJson = prettyjson(jsonencode(init));
                        fprintf(fid,'%s',formattedJson);
                        fclose(fid);
                    end
                otherwise
                    error(['Output type ''',type,''' not supported'])
            end
        end

    end
end