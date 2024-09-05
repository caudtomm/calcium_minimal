function registration_opticflow_singlestack(FilenameRAW,ref_img)


FilenameRAW = char(FilenameRAW);

fspecs = getFileNameSpecs(FilenameRAW);

disp('')
disp('OpticFlow registration:')
disp(fspecs.fname)

% output directory
FileOut_path = fullfile(fspecs.orig_fpath, 'flow_registration_results');

%%
options = OF_options(...
    'input_file', FilenameRAW, ...
    'output_path', fullfile(FileOut_path,fspecs.fname), ... % results folder
    'output_format', 'TIFF', ... % output file format: HDF5, MAT, 
                             ... % TIFF or MULTIFILE_HDF5, ... to generate multiple files
                             ... % or CAIMAN_HDF5 for CAIMAN support
    'alpha', 10, ... % smoothness parameter
    'sigma', [1, 0.5, 0.1; ...  % gauss kernel size channel 1
              1, 0.5, 0.1], ... % gauss kernel size channel 2
    'quality_setting', 'quality', ... % set the quality out of 'fast', 'medium' or 'quality' (default)
    'bin_size', 1, ... % binning over 1 frames from the data
    'buffer_size', 500, ... % size of blocks for the parallel evaluation (larger takes more memory)
    'reference_frames', ref_img ...
    );

% saving the options to txt file (for archiving):
options.save_options(fullfile(FileOut_path, 'options.json'));


compensate_recording(options);

cd(FileOut_path)
end
