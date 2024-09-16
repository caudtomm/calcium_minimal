classdef OpticFlowRegistration < Registration
    properties
        % Re-instantiate inherited properties
        data_processed Movie
    end

    methods (Access=protected)
        function [obj,movie_result,operation] = registration_alg(obj, data_raw)
            arguments
                obj OpticFlowRegistration
                data_raw Movie
            end

            % initialize
            movie_result = data_raw;
            movie = data_raw.stack;
            ref_img = obj.reference_img;
            if isempty(ref_img); ref_img = mean(movie,3,'omitmissing'); end
            obj.reference_img = ref_img;
            

            options = OF_options(...
                'input_file', FileIn, ...
                'output_path', FileOut, ... % results folder
                'output_format', 'TIFF', ... % output file format: HDF5, MAT, 
                                         ... % TIFF or MULTIFILE_HDF5, ... to generate multiple files
                                         ... % or CAIMAN_HDF5 for CAIMAN support
                'alpha', 10, ... % smoothness parameter
                'sigma', [1, 0.5, 0.1; ...  % gauss kernel size channel 1
                          1, 0.5, 0.1], ... % gauss kernel size channel 2
                'quality_setting', 'quality', ... % set the quality out of 'fast', 'medium' or 'quality' (default)
                'bin_size', 1, ... % binning over 1 frames from the data
                'buffer_size', 500, ... % size of blocks for the parallel evaluation (larger takes more memory)
                'reference_frames', obj.reference_img ...
                );
            
            % saving the options (for archiving):
            obj.init.options = options;
                        
            compensate_recording(options);
        end

    end
    
    methods
        % Constructor
        function obj = OpticFlowRegistration(movie_in, reference_img)
            arguments
                movie_in = ''
                reference_img double = []
            end
            obj = obj@Registration(movie_in,reference_img);

            % definition of init structure
            obj.init.method = 'OpticFlow Registration';

            % initialization
            obj.data_processed = obj.data_raw;
            obj.operation = [];
        end

    end
end