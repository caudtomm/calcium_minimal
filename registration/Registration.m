classdef (Abstract) Registration < MovieProcessing
    % My description
    properties
        operation
        reference_img double
    end
    
    methods (Abstract, Access=protected)
        registration_alg
    end

    methods
        % Constructor
        function obj = Registration(movie_in,reference_img,reference_meta)
            arguments
                movie_in = ''
                reference_img double = []
                reference_meta struct = struct()
            end
            obj = obj@MovieProcessing(movie_in);
            obj.reference_img = reference_img;
            obj.init.reference_meta = reference_meta;
            
        end

        %% Funtional methods
        function obj = run(obj)
            % Signal registration start
            obj.disp_runheader()

            % initialize
            movie = obj.data_raw;

            % Remove bad periods
            rbp = BasicMovieProcessor('remove_badperiods',movie);
            rbp = rbp.run;
            movie = rbp.data_processed;
            
            % register
            [obj, movie, op] = obj.registration_alg(movie);
            obj.data_processed = movie;
            obj = obj.tailsequence();

            % add nans in the place of bad periods
            abp = BasicMovieProcessor('add_badperiods_as_nans',movie);
            abp = abp.run;
            movie = abp.data_processed;
            
            % store results
            movie = movie.checkAgain;
            obj.data_processed = movie;
            obj.operation = op;
        end


    end

end