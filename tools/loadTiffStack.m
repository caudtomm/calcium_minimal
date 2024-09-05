function [stack, headers] = loadTiffStack(tifffilename,varargin)


p = inputParser;
if verLessThan('matlab','8.2')
    p.addParamValue('firstIdx',1,@isnumeric)
    p.addParamValue('lastIdx',1,@isnumeric)
    p.addParamValue('option',1,@isnumeric)
    p.parse(varargin{:})
else
    addRequired(p,'tifffilename',@ischar)
    addParameter(p,'firstIdx',1,@isnumeric)
    addParameter(p,'lastIdx',1,@isnumeric)
    addParameter(p,'option',1,@isnumeric)
    parse(p,tifffilename,varargin{:})
    tifffilename = p.Results.tifffilename;
end

if nargout > 1
    loadHeaders = true;
else
    loadHeaders = false;
end

fprintf('Loading TIFF frame ...');

warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');
warningsBackOn = onCleanup(...
    @() warning('on', 'MATLAB:imagesci:tiffmexutils:libtiffWarning'));

if ischar(tifffilename)
    tiff = Tiff(tifffilename, 'r');
    closeTiff = onCleanup(@() close(tiff));
end

if p.Results.option == 1
    w = tiff.getTag('ImageWidth');
    h = tiff.getTag('ImageLength');
    dataClass = class(read(tiff));
    nFrames = findNoOfFrames(tiff,1001);
    stack = zeros(h, w, nFrames, dataClass);
    if loadHeaders
        headers = cell(1, nFrames);
    end
    setDirectory(tiff, p.Results.firstIdx);
    for t = 1:nFrames
        stack(:,:,t) = read(tiff);
        
        if loadHeaders
            headerNames = tiff.getTagNames;
            headers{t} = getTag(tiff, 'ImageDescription');
            try
                headers{t} = [headers{t} getTag(tiff,'Software')];
            catch
            end
        end
        
        if t < nFrames
            nextDirectory(tiff);
        end
    end
    
elseif p.Results.option == 2
    info = imfinfo(tifffilename);
    offset = info(1).Offset;
    w = info(1).Width;
    h = info(1).Height;
    dataClass = 'int16'; % this is true for our ScanImage recordings, 
    % it is possible to use info.BitDepth to try and figure out nBytesPerSample
    
    % MK: using memmapfile here, but just reading as a binary file might be
    % faster, worth trying
    m = memmapfile(tifffilename, 'Format', dataClass, 'Offset', offset);
    data = m.Data;
%     clear m;
    nPixels = w*h;
    frameIdx = p.Results.firstIdx:1:length(info);
    data = reshape(data, [], length(info));
    nSamples = size(data, 1); % number of values in each frame related vector
    % the images are the last nPixels values of this vector
    % here we assume all the headers occupy the same number of bytes,
    % this is true for ScanImage recordings. info.Offset is a useful field
    % otherwise.
    stack = data(nSamples-nPixels+1:nSamples, frameIdx);
    stack = reshape(stack, w, h, []);
    stack = permute(stack, [2 1 3]);
    if loadHeaders
        headers = {info(frameIdx).ImageDescription};
    end
end

fprintf(' done\n');

end % of main function


function n = findNoOfFrames(tiff, seekInterval)
if nargin < 2
    seekInterval = 1000;
end
%keep guessing until we seek too far
guess = seekInterval;
overSeeked = false;
n = 0;

while ~overSeeked
    try
        tiff.setDirectory(guess);
        guess = 2*guess; %double the guess
    catch ex
        overSeeked = true; %we tried to seek past the last directory
    end
end
%when overseeking occurs, the current directory/frame will be the last one
n = tiff.currentDirectory;
end


