classdef (Abstract) SingleProcess < Process
    properties (Abstract)
        operation
    end

    methods (Abstract, Access=protected)
        tailsequence % sequence of events upon Process run() completion
    end

    methods (Access=protected)
        function disp_runheader(obj)
            disp('')
            disp('')
            disp(['Run started: ', obj.init.method]);
            if ~isempty(obj.data_raw) && isprop(obj.data_raw, 'path') && ~isempty(obj.data_raw.path)
                disp(obj.data_raw.path.fname);
            end
            disp('')
        end       
    end
    
    methods
        function obj = SingleProcess()

        end

    end
end