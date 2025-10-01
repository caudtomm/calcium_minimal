function outfiles = cutMovie (folder,varargin)
% cutMovie (folder,varargin)
% Cuts movies of a certain format (default = .avi) using ffmpeg in the folder
% specified. When multiple cutouts are requested, they are placed back-to-back.
%
% Optional: 'start' (sec), 'durSec' (sec), 'inForm', 'outForm', 'outDir', 'nOut'
% Example:
% cutMovie(f, 'start', 0, 'durSec', 3600, 'inForm', '.avi', 'outForm', '.avi', 'nOut', 10)

p = inputParser;
p.KeepUnmatched = true;

defaultinForm = '.avi';
defaultDurSec = 7200; % 2h
defaultStartT = 0;
defaultOutdir = pwd;

addRequired(p,'folder',@ischar);
addParamValue(p,'inForm',defaultinForm,@ischar)
addParamValue(p,'start',defaultStartT,@isnumeric)
addParamValue(p,'durSec',defaultDurSec,@isnumeric)
addParamValue(p,'outDir',defaultOutdir,@ischar)

parse(p,folder,varargin{:})

durSec = p.Results.durSec;
outdir = p.Results.outDir;

cd(folder)
files = dir(['*',p.Results.inForm]);
outfiles = cell(numel(files),1);

% Common ffmpeg quiet flags
ffquiet = '-hide_banner -loglevel error -nostdin -y';

for kk = 1:length(files)
    try
        info = mmfileinfo(files(kk).name);
    catch
        info  = VideoReader(files(kk).name);
    end
    infiledur = info.Duration; % sec

    if isfield(p.Unmatched, 'dur')
        durSec = p.Unmatched.dur;
    end

    if isfield(p.Unmatched, 'nOut')
        nc = p.Unmatched.nOut;
        if nc > ceil(infiledur/durSec) && ~isfield(p.Unmatched, 'dur')
            durSec = round(infiledur/nc);
        elseif isfield(p.Unmatched, 'dur') && nc > ceil(infiledur/p.Unmatched.dur)
            nc = ceil(infiledur/durSec);
        end
    else
        nc = ceil(infiledur/durSec);
    end

    duration = datestr(durSec/(3600*24), 'HH:MM:SS');

    if isfield(p.Unmatched, 'outForm')
        outForm = p.Unmatched.outForm;
    else
        outForm = p.Results.inForm;
    end

    for mm = 1:nc
        stSec = (p.Results.start + (mm-1)*durSec);
        st = datestr(stSec/(3600*24),'HH:MM:SS');

        endstr  = ['_p',sprintf('%02d',mm),outForm];
        outname = strrep(files(kk).name,p.Results.inForm,endstr);
        outname = fullfile(outdir,outname);
        outfiles{kk} = outname;

        % Minimal progress line
        fprintf('[%d/%d] %s -> %s  (start %s, dur %s)\n', ...
            mm, nc, files(kk).name, outname, st, duration);

        inq  = ['"', files(kk).name, '"'];
        outq = ['"', outname, '"'];

        if ~strcmp(p.Results.inForm,outForm)
            cmd = sprintf('ffmpeg %s -i %s -ss %s -t %s %s', ...
                ffquiet, inq, st, duration, outq);
        else
            cmd = sprintf('ffmpeg %s -i %s -ss %s -t %s -vcodec copy %s', ...
                ffquiet, inq, st, duration, outq);
        end

        % Run quietly; only show something if it fails
        [status,cmdout] = system(cmd);
        if status ~= 0
            warning('ffmpeg failed: %s\nCommand:\n%s', strtrim(cmdout), cmd);
        end
    end
end
