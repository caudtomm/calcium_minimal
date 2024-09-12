classdef (Abstract) MovieProcessing < SingleProcess
    % Any operation performed on a Tiff stack (2p movie)
    properties
        data_raw Movie = Movie('') % setting directly is discouraged - if possible use method 'setDataRaw()' for more robust performance 
    end

    properties (Abstract)
        data_processed Movie % specifying accepted data type
    end

    methods (Access=protected)
        function obj = tailsequence(obj,log)
            arguments
                obj
                log char = class(obj)
            end

            movie = obj.data_processed;

            disp('')
            movie = movie.checkAgain;
            newlog = [obj.hash,' : ',log];
            movie = movie.update_log(newlog);
            obj.data_processed = movie;
        end
    end

    methods
        function obj = MovieProcessing(movie_in)
            if exist('movie_in','var') && ~isempty(movie_in)
                if isa(movie_in, 'Movie')
                    obj.data_raw = movie_in;
                elseif ischar(movie_in) || isstring(movie_in)
                    % Assume movie_in is a path to a .mat file containing a 'Movie' object
                    obj.data_raw = load(movie_in, 'Movie'); % Ensure 'Movie' is the variable name in the .mat file
                else
                    error('Input must be of type Movie or a valid path to a .mat file containing a Movie object.');
                end
            end
        end
        
        function obj = setDataRaw(obj, value)
            if isa(value,'Movie')
                obj.data_raw = value;
            else
                error('Data must be a Movie object')
            end
        end

    end
end