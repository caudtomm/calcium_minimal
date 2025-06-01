classdef ComponentReconstruction < MovieProcessing
    properties
        data_processed Movie = Movie('')
        operation
    end

    methods
        function obj = ComponentReconstruction(operation, movie_in)
            arguments
                operation 
                movie_in = ''
            end
            obj = obj@MovieProcessing(movie_in);
            obj.operation = operation;
        end

        function obj = run(obj)
            % initialize
            inputdata = obj.data_raw;
            outputdata = inputdata; % init for badperiods
            pxmu = obj.operation.pxmu;
            pxstd = obj.operation.pxstd;
            W = obj.operation.transform;
            idx_nanpx = obj.operation.idx_nanpx;
            method = obj.operation.inputmethod;
            mode = obj.operation.mode;
            if isempty(mode); mode = 'nooffset'; end
            do_subtr = obj.operation.subtract_from_inputdata;
            
            fprintf('Reconstructing movie component - mode: %s\n',mode)

            if do_subtr
                % save a copy of the original data;
                orig_data = inputdata.stack;
            end

            % linearize input data
            inputdata = BasicMovieProcessor('remove_badperiods',inputdata).run().data_processed;
            disp('Linearizing frames ...')
            inputdata = PCAICA.linearizeFrames(inputdata.stack,true,[],idx_nanpx,method);
            
            % reconstruct selected IC
            disp('Reconstructing component ...')
            switch mode
                case 'nooffset'
                    reconstructed = (W * ((inputdata-pxmu)./pxstd)) .* pxstd;
                case 'withoffset'
                    reconstructed = (W * ((inputdata-pxmu)./pxstd)) .* pxstd + pxmu;
                otherwise
                    error('specified reconstruction mode unknown.')
            end

            % check how to properly deal with idx_nanfr

            % recontruct full movie
            movie_fullyrc = PCAICA.addnanpxsbackin(reconstructed, idx_nanpx);
            movie_fullyrc = reshape(movie_fullyrc,obj.data_raw.h,obj.data_raw.w,[]);

            % add badperiods back in
            outputdata.stack = movie_fullyrc;
            outputdata = BasicMovieProcessor('add_badperiods_as_nans',outputdata).run().data_processed;

            % get some space back for subtration, but also you don't wanna
            % save a massive operation into the init mat file
            obj.operation = [];
            clear W

            % subtract from input data, if requested
            if do_subtr
                disp('Subtracting from input movie')
                outputdata.stack = orig_data - outputdata.stack;
            end

            % store to results
            obj.data_processed = outputdata;
        end
    end
end