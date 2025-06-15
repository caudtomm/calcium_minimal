function save_avi_crop(inputFile, outputFile, cropRect, frameRange)
% save_avi_crop(inputFile, outputFile, cropRect, frameRange)
% Crops a region from a large AVI file and saves it without loading full video into memory.
% 
% Parameters:
%   inputFile  - path to input AVI file
%   outputFile - path to save cropped AVI
%   cropRect   - [x, y, width, height] for cropping (pixels)
%   (optional) frameRange - [startFrame, endFrame] to crop (1-based indexing)

v = VideoReader(inputFile);

% Validate frame range
if nargin>3
    frameRange(1) = max(1, frameRange(1));
    frameRange(2) = min(v.NumFrames, frameRange(2));
else
    frameRange = [1 v.NumFrames];
end

% create output folder if necessary
outfolder = getFileNameSpecs(outputFile).orig_fpath;
if ~exist(outfolder,"dir"); mkdir(outfolder);end

% Setup VideoWriter
vw = VideoWriter(outputFile, 'Grayscale AVI');
vw.FrameRate = v.FrameRate;
open(vw);

frameCount = 0;

for k = frameRange(1):frameRange(2)
    v.CurrentTime = (k-1) / v.FrameRate;
    if hasFrame(v)
        frame = readFrame(v);
        % Crop using imcrop or indexing
        x = cropRect(1); y = cropRect(2); w = cropRect(3); h = cropRect(4);
        crop = frame(y:y+h-1, x:x+w-1, :);
        writeVideo(vw, rgb2gray(crop));
        frameCount = frameCount + 1;
    else
        warning('End of video reached at frame %d.', k);
        break;
    end
end

close(vw);
fprintf('Saved %d cropped frames to %s\n', frameCount, outputFile);
end
