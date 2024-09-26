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

        function [v,h,collage] = subsampledMovieCollage(obj)
            % retrieve file list
            fileList = {};
            for i = 1:numel(obj.sj.filelist)
                fileList{end+1} = find_daughter_file(...
                    fullfile(obj.folder,obj.sj.filelist(i).name),'mat');
            end

            % Determine the number of files
            numFiles = length(fileList);
        
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
            
            % Create a figure to display the collage
            h = figure;
            colormap(gray);
        
            v = VideoWriter(fullfile(obj.folder,'sidebyside.avi'));
            open(v);

            maxval = quantile(collage(:),.999);

            % Play the collage frame by frame
            for frame = 1:z
                imagesc(collage(:, :, frame),[0 maxval]);
                axis image;
                title(['frame: ',num2str(frame)])
                normframe = collage(:, :, frame)/maxval;
                normframe(normframe>1) = 1;
                normframe(normframe<0) = 0;
                writeVideo(v,normframe);
                % pause(1/16); % 16Hz playback
            end

        end

    end
end
