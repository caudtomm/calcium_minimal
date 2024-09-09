classdef OpticFlowRegistration < Registration
    properties
        % Re-instantiate inherited properties
        data_processed Movie
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
            
            % saving the options to txt file (for archiving):
%             options.save_options(fullfile(FileOut_path, 'options.json'));
            obj.init.options = options;
                        
            compensate_recording(options);
            
            
%             if isempty(shifts)
%                 disp('Computing shifts..')
%                 [shift_x, shift_y] = register_frames(movie, ref_img);
% %                 shift_x = mean(shift_x) * ones(size(shift_x));
% %                 shift_y = mean(shift_y) * ones(size(shift_y));
%             else
%                 disp('Applying input shifts..')
%                 shift_x = shifts(:,1);
%                 shift_y = shifts(:,2);
%             end
%             movie = shift_data(movie,shift_x,shift_y);
%             
%             movie_result.stack = movie;
%             operation.shifts = [squeeze(shift_x),squeeze(shift_y)];
        end
    end
end