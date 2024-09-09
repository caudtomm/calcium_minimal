% Peter Rupprecht 02-2015 using Scanimage B data
% reads metadata and returns image info file and framerate, zstep,
% framerate, zoom, motorpositions, scalingfactors
% optional input is the -full- filename/path
function [A,result,framerate,zstep,zoom,motorpositions,scalingfactors] = read_metadata_function(filename)

if strcmp(filename,'')
    %% read in path/file/folder of image
    PathName = '';
    [FileName,PathName,FilterIndex] = uigetfile(strcat(PathName,'*.tif')); file_name = strcat(PathName,FileName);
else
    file_name = filename;
end
A = imfinfo(file_name);


%% define metadata of interest 

% snippet{1} = 'state.acq.frameRate';
% snippet{2} = 'state.internal.triggerFrameDelayMS';
% snippet{3} = 'state.internal.triggerTimeFirstString';
% snippet{4} = 'state.acq.zStepSize';
% snippet{5} = 'state.acq.numberOfZSlices';
% snippet{6} = 'state.acq.pixelsPerLine';
% snippet{7} = 'state.acq.zoomFactor';
% snippet{8} = 'state.acq.linesPerFrame';
% snippet{9} = 'state.acq.inputBitDepth';
% snippet{10} = 'state.internal.triggerTimeString';
% snippet{11} = 'state.acq.frameRate';
% snippet{12} = 'state.acq.frameRate';
% snippet{13} = 'state.init.eom.maxPower';
% snippet{14} = 'state.acq.numAvgFramesSaveGUI';
% snippet{15} = 'state.acq.numAvgFramesSaveGUI';
% snippet{16} = 'state.acq.numberOfFrames';
% snippet{17} = 'state.acq.numAvgFramesSaveGUI';

snippet{1} = 'SI.hRoiManager.scanFrameRate';
snippet{2} = 'state.internal.triggerFrameDelayMS';
snippet{3} = 'state.internal.triggerTimeFirstString';
snippet{4} = 'SI.hStackManager.stackZStepSize';
snippet{5} = 'SI.hStackManager.numSlices';
snippet{6} = 'SI.hRoiManager.pixelsPerLine';
snippet{7} = 'SI.hRoiManager.scanZoomFactor';
snippet{8} = 'SI.hRoiManager.linesPerFrame';
snippet{9} = 'state.acq.inputBitDepth';
snippet{10} = 'state.internal.triggerTimeString';
snippet{11} = 'SI.hRoiManager.scanFrameRate';
snippet{12} = 'SI.hRoiManager.scanFrameRate';
snippet{13} = 'state.init.eom.maxPower';
snippet{14} = 'state.acq.numAvgFramesSaveGUI';
snippet{15} = 'state.acq.numAvgFramesSaveGUI';
snippet{16} = 'SI.hStackManager.framesPerSlice';
snippet{17} = 'state.acq.numAvgFramesSaveGUI';

snippet{18} = 'state.motor.absZPosition';
snippet{19} = 'state.init.scanOffsetAngleX';
snippet{20} = 'SI.hRoiManager.scanFrameRate';

%% read out metadata
clear result
% reference = A(1).ImageDescription;
reference = A(1).Software;
% lag = 1;
lag = 3;
for jj = 1:17
    k = strfind(reference(1:end),snippet{jj});
    check = 0;
    stringlength = 10;
    while check == 0
    result{jj} = str2num( reference(k + length(snippet{jj})+lag:k + length(snippet{jj}) + lag + stringlength) );
    if isempty(result{jj}) && stringlength>0
        stringlength = stringlength - 1;
    else
        if stringlength<=0; warning(strcat('Snippet #',num2str(jj),' : not found')); end
        check = 1;
    end
    end
end

% exceptions

k = strfind( reference(1:end),snippet{18});
result{18} = reference(k + length(snippet{18})+lag:k + length(snippet{18}) + lag + 24);
k = strfind( reference(1:end),snippet{19});
try
    result{19} = reference(k + length(snippet{19})+lag:k + length(snippet{19}) + lag + 19);
catch
    result{19} = reference(k + length(snippet{19})+lag:end);
end
k = strfind( reference(1:end),snippet{20});
try
    result{20} = reference(k + length(snippet{20})+lag:k + length(snippet{20}) + lag + 5);
catch
    result{20} = [];
end
%% make metadata it easier to process

if result{12}; framerate = min(result{11},result{1}); else framerate = result{1}; end
if ~isempty(result{20}); framerate = str2double(result{20}); end
zstep = result{4};
zoom = result{7};
motorpositions = result{18};
scalingfactors = result{19};

end