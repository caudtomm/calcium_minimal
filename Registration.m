classdef Registration < MovieProcessing
    % My description
    properties
        method
        reference_img
        init
    end
    
    methods
        % Constructor
        function obj = Registration(filename)
%             obj.movie_raw = double(loadTiffStack(filename));
% 
%             obj.init.filename = filename;
%             obj.init.meta = get_meta(filename);

        end

        %% Getters / setters

        function meta = get_meta(Filename)
            L = imfinfo(Filename);
            meta.height = L(1).Height;
            meta.width = L(1).Width;
   
            meta.numberframes = size(movie,3);
        end

        %% Funtional methods

        function run()
            %% Signal registration start
            disp_runheader()

            %% Subtract bad periods
            [movie, ~] = subtract_badperiods(movie,badperiods);
            
            %% register
            registration_alg()
%             [movie,movie_raw,shift_x,shift_y] = linearAlignMovie(movie,movie_raw,ref_img);
            
            %% add nans in the place of bad periods
            movie = add_badperiods_as_nans(movie,badperiods,meta);
        end

        function registration_alg()
            error('Default behavior unspecified.')
        end

        function disp_runheader()
            disp('')
            disp(['Run started: ',method,' registration:'])
            disp(fspecs.fname)
        end

    end
end