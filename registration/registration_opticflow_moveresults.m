function registration_opticflow_moveresults(loc,do_correct)

% go to results folder
cd(loc.warp.rawtrials_opticflow_unchecked);

% create output folder

outputfolder = loc.warp.rawtrials_opticflow;
FileOut_path = fullfiletol(loc.subject_datapath,outputfolder);

mkdir(FileOut_path);

%% move all optic-flow-compensated stacks to target directory

filelist = dir();
allbadperiods = [];

for i_f = 1:numel(filelist)
    filename = filelist(i_f).name;

    % skip files and parent directories
    if ~filelist(i_f).isdir; continue; end
    if ismember(filename,{'.','..'}); continue; end

    % manage extensions
    foldername = filename;
    FileIn_ext = loc.datafile_ext;
    if endsWith(foldername,FileIn_ext)          % rename folder to exclude extension
        foldername_new = foldername(1:length(foldername)-length(FileIn_ext)-1);
        movefile(foldername,foldername_new)
        foldername = foldername_new;
    else                                        % make sure file name includes extension
        filename = strcat(filename,'.',FileIn_ext);
    end
    
    % prepare input and output names
    statsname = fullfiletol(foldername,'statistics.mat');
    inputname = fullfiletol(foldername,'compensated.tiff');
    outputname = fullfiletol(FileOut_path,filename);
    

    disp(filename)

    if do_correct
        % load registration statistics, identify bad periods and replace
        % them with a different source: typically with the rigid
        % registration results.
        thisbw = replace_badwarp_opticflow(loc,inputname,outputname,statsname);
        allbadperiods = [allbadperiods; [repelem(i_f,size(thisbw,1),1), thisbw]];
    else
        % copy file
        copyfile(inputname,outputname)
    end

end

saveMatrixToCSV(allbadperiods, ...
    fullfiletol(outputfolder,'auto_badperiods.csv'), ...
    {'relative trial num', ...
    'bad period start frame (included)', ...
    'bad period end frame (included)'});

end



function thisbw = replace_badwarp_opticflow(loc,inputname,outputname,statsname)
    loc.warptag = '';
    loc.fspecs = getFileNameSpecs(inputname);
    loc.rawpath = fullfiletol(loc.subject_datapath,loc.rawtrials_rigidreg);
    loc.WarpConnectPath = fullfiletol( ...
        loc.rawpath,loc.warp.rawtrials_opticflow_unchecked, ...
        loc.fspecs.orig_fpath);
    loc.rawfname = loc.fspecs.orig_fpath;

    %%Change this according to your imagej
    javaaddpath([loc.drive,'\scratch\gfriedri\caudtomm\code\Fiji.app\jars\ij-1.53t.jar']);
        
    % select bad periods automatically
    thisbw = autoselectBadWarp_opticflow(statsname);
    if isempty(thisbw); copyfile(inputname,outputname); return; end

    % replace identified bad period with the corresponding periods in the
    % 'rawpath' (normally, rigidly registered raw stacks - or anyways,
    % whatever was fed into the warping algorhythm in the first place 
    replaceBadWarp(inputname,outputname,thisbw,loc)

end


function thisbw = autoselectBadWarp_opticflow(statsname)
    % load the chosen metric from statistics file
    mymetric = load(statsname).mean_translation;

    % strategy: during bad periods, the mean translation is stuck at a
    % static value. Identify this feature by taking both the fist and
    % second derivative = 0 (to ensure that it is not a single value=0).
    % The first derivative's 'i' value corresponds to 'i+1' in the original
    % vector; but the first 'equal' value in the original vector
    % corresponds to the prior value to the first i=0 in the first
    % derivative. These two effects cancel out, and we can take the indices
    % where dS/dt = 0 directly as indices of the badly warped frames.
    idx = diff(mymetric,1)==0 & [0,diff(mymetric,2)]==0;

    % identify ranges
    thisbw = findContinuousStretches(idx);

end


function saveMatrixToCSV(matrix, filename, header)
% saveMatrixToCSV - Save a matrix to a CSV file with an optional header.
%
% Syntax:
%   saveMatrixToCSV(outputMatrix, filename, header)
%
% Description:
%   This function saves a numeric matrix to a CSV (Comma-Separated Values)
%   file. You can provide an optional header to include column names in the
%   CSV file.
%
% Inputs:
%   outputMatrix - The numeric matrix to be saved to the CSV file.
%   filename     - The name of the CSV file to create or overwrite.
%   header       - (Optional) A cell array of column names for the CSV
%                  file. Default is an empty cell array.
%
% Example:
%   outputMatrix = [1, 5; 7, 10; 15, 20];
%   filename = 'output.csv';
%   header = {'Start Index', 'End Index'};
%   saveMatrixToCSV(outputMatrix, filename, header);
%
% See also:
%   -
%
% Author: Tommaso Caudullo
% Date: 17/09/2023

    % Check the number of input arguments
    narginchk(2, 3);
    
    % Default header if not provided
    if nargin < 3
        header = {};
    end

    % Add the header to the output matrix
    data = [header; num2cell(matrix)];

    % Write the data to the CSV file
    try
        fid = fopen(filename, 'w');
        if fid == -1
            error('Unable to create the CSV file.');
        end

        % Write each row to the file
        for i = 1:size(data, 1)
            fprintf(fid, '%s,', data{i, 1:end-1});
            fprintf(fid, '%s\n', data{i, end});
        end

        % Close the file
        fclose(fid);
        disp(['CSV file saved: ', filename]);
    catch
        if exist('fid', 'var')
            fclose(fid);
        end
        error('Error writing to the CSV file.');
    end
end


