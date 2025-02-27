classdef RegistrationViewer
    properties
        sj
        folder
    end

    methods
        function obj = RegistrationViewer(subject, foldername)
            arguments
                subject Subject
                foldername char = pwd
            end

            obj.sj = subject;            
            obj.folder = foldername;
        end

        function collage = subsampledMovieCollage(obj)
            orig_folder = fullfiletol(pwd);

            % get file list and move to subject-datapath
            [fileList, numFiles] = obj.sj.getFileListInSubfolder(obj.folder);
        
            % Calculate the dimensions of the collage
            cols = ceil(sqrt(numFiles));
            rows = ceil(numFiles / cols);
        
            % Load the first file to get the size of movie.stack
            disp('Initializing...')
            load(fileList{1}, 'movie');
            [x, y, z] = size(movie.stack);
            
            % Subsample dimensions
            x = ceil(x / 5);
            y = ceil(y / 5);
        
            % Initialize the collage matrix
            collage = zeros(rows * x, cols * y, z);
            
            % Add each movie to the collage
            for i = 1:numFiles
                % Load the movie data
                disp(fileList{i})
                load(fileList{i}, 'movie');
                
                % Subsample the movie stack
                subsampledStack = movie.stack(1:5:end, 1:5:end, :);
                
                % Determine the position in the collage
                rowIdx = ceil(i / cols);
                colIdx = mod(i-1, cols) + 1;
                
                % Insert the subsampled movie stack into the collage
                collage((rowIdx-1)*x + (1:x), (colIdx-1)*y + (1:y), :) = subsampledStack;
            end

            collage = Movie(collage);
            collage = collage.setFramerateHz(obj.sj.framerate);
            collage.path.fname = 'sidebyside';
            
            % return to original location
            cd(orig_folder)
        end

        function [fbf_corr, zsc_corr, putative_fails] = QC(obj)
            % generate a frame-by-frame correlation curve

            orig_folder = fullfiletol(pwd);

            % get file list and move to subject-datapath
            [fileList, numFiles] = obj.sj.getFileListInSubfolder(obj.folder);
            if numFiles~=obj.sj.getNTrials; warning('Number of movies found is different from the number of trials.'); end

            % initialize
            fbf_corr = zeros(obj.sj.getNFrames-1, numFiles); % [numFrames-1, trials]

            % loop
            parfor i = 1:numFiles
                % Load the movie data
                disp(fileList{i})
                m = load(fileList{i}, 'movie');

                fbf_corr(:,i) = avg_frame_correlation(m.movie.stack);
            end

            % zscored frame-by-frame correlation
            zsc_corr = nanzscore(fbf_corr);

            % define putatively failed frames
            putative_fails = [zeros(1,numFiles) ; zsc_corr < -5]; % add one line at the start for frame 1

            % return to original location
            cd(orig_folder)
        end

    end
end


function corr_curve = avg_frame_correlation(data)
    [rows, cols, frames] = size(data);
    corr_curve = zeros(frames - 1, 1);
    
    for i = 1:frames-1
        frame1 = reshape(data(:,:,i), [], 1);
        frame2 = reshape(data(:,:,i+1), [], 1);
        corr_curve(i) = corr(frame1, frame2);
    end
end