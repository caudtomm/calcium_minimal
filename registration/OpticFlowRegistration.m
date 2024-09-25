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

            % define reference
            movie = data_raw.stack;
            ref_img = obj.reference_img;
            if isempty(ref_img); ref_img = mean(movie,3,'omitmissing'); end
            obj.reference_img = ref_img;

            % create a temporary tiff file (typically useful to have no nans)
            data_raw.path.fname = ['temp',data_raw.path.fname];
            TiffFilename = data_raw.save(obj.init.original_path,'tif');

            % name the output tiff file
            fpathout = obj.hash;
            if obj.init.batch_run; fpathout = ['batch_',obj.init.superHash]; end
            fpathout = fullfile(obj.init.original_path,fpathout,movie_result.path.fname);

            % flow registration settings
            options = OF_options(...
                'input_file', TiffFilename, ...
                'output_path', fpathout, ... % results folder
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
            
            % saving the settings (for archiving):
            obj.init.internal_settings = options;
                        
            % run flow registration
            compensate_recording(options);

            % load result stack into movie
            movie_result.stack = Movie(fullfile(fpathout,'compensated.TIFF')).stack;

            % load operation results
            operation = load(fullfile(fpathout,'statistics.mat'));

            % clean up: delete temporary tiff file
            delete(TiffFilename)
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