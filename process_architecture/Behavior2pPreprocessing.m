classdef Behavior2pPreprocessing
    % Class for preprocessing 2-photon behavior AVI files with ffmpeg.
    % MJPEG encoding.
    
    properties
        ffm = 'ffmpeg'        % path to ffmpeg executable
        outdir = ''           % default output directory ('' = same as input)
        suf = ''              % default suffix for rotated files
        q = 3                 % ffmpeg MJPEG quality parameter
        pix = 'yuvj422p'      % pixel format for ffmpeg
    end  
        
    methods (Static)
        function out = cut_avis_ffmpeg(src, dur_s, outdir)
            % Network connection is the bottleneck
            out = cutMovie(src,'durSec',dur_s, 'outDir', outdir);
            cd(outdir)
        end

        function cmd = run_fiji_batch(subjectID, rois, ijExe, ijmFile, workers, rankRange)
        %RUN_FIJI_BATCH Launches Fiji macro over multiple ranks with ROI args from a Map.
        %
        % Inputs
        %   subjectID  : char/string, passed as subjectID=... to the macro
        %   rois       : containers.Map with values having fields Xs (1x2), Ys (1x2)
        %                e.g. rois('head').Xs = [x1 x2]; rois('head').Ys = [y1 y2]
        %   ijExe      : (optional) full path to ImageJ-win64.exe
        %   ijmFile    : (optional) full path to the .ijm macro file
        %   workers    : (optional) integer workers=...
        %   rankRange  : (optional) [start step end] for FOR /L (default [0 1 5])
        %
        % Output
        %   cmd        : the full command string that was executed (for logging)
        %
        % Notes
        %   - Each ROI becomes KEY=x,y,w,h with KEY = upper(name), ordered by rois.keys.
        %   - x,y,w,h are computed from [min(Xs) min(Ys) diff(Xs)+1 diff(Ys)+1].
        %   - Uses START to spawn a process per rank.
        %
        % Example
        %   cmd = run_fiji_batch('TC_230506_TC0003_230502beh2b2_sxpDp_odorexp004_RPB3144501500AG', ...
        %                        rois, ...
        %                        'W:\scratch\gfriedri\caudtomm\code\Fiji.app\ImageJ-win64.exe', ...
        %                        'W:\scratch\gfriedri\caudtomm\code\Fiji.app\macros\random_IJmacro_scripts\batch_tail_movies4.ijm', ...
        %                        6, [0 1 5]);
        
            if nargin < 3 || isempty(ijExe)
                ijExe = 'W:\scratch\gfriedri\caudtomm\code\Fiji.app\ImageJ-win64.exe';
            end
            if nargin < 4 || isempty(ijmFile)
                % this program cannot run headless (Log saving)
                ijmFile = 'W:\scratch\gfriedri\caudtomm\code\Fiji.app\macros\random_IJmacro_scripts\batch_tail_movies5.ijm';
            end
            if nargin < 5 || isempty(workers)
                workers = 6;
            end
            if nargin < 6 || isempty(rankRange)
                rankRange = [0 1 5];
            end
        
            ks = rois.keys;
            ks = sort(lower(string(ks)));                             % stable order, case-insensitive
            roiArgs = strings(1, numel(ks));
            for i = 1:numel(ks)
                k = ks(i);
                r = rois(char(k));
                x0 = min(r.Xs(:));
                y0 = min(r.Ys(:));
                w  = abs(diff(r.Xs(:)))./1; w = w(1) + 1;             % width in pixels
                h  = abs(diff(r.Ys(:)))./1; h = h(1) + 1;             % height in pixels
                rect = sprintf('%d,%d,%d,%d', round(x0), round(y0), round(w), round(h));
                roiArgs(i) = upper(k) + "=" + rect;                   % e.g., HEAD=249,486,56,44
            end
        
            macroArgs = strjoin( ...
                ["workers=" + string(workers), ...
                 "rank=%R", ...
                 "subjectID=" + string(subjectID), ...
                 roiArgs], " ");
        
            % Build the Windows FOR /L command; avoid sprintf to keep '%' literal.
            forSpec   = "(" + strjoin(string(rankRange), ",") + ")";
            titleStr  = '""';  % empty START title
            
            cmd = "for /l %R in " + forSpec + " do " + ...
                  "start " + titleStr + " " + quote(ijExe) + ...
                  " --ij2 -macro " + quote(ijmFile) + " " + quote(macroArgs);
        
            % Execute; capture status if you want (status==0 means success).
            [status, out] = system(cmd); %#ok<ASGLU> 
            if status ~= 0
                warning('run_fiji_batch:systemFailed', 'Non-zero exit from system(). Inspect your paths and arguments.');
            end
        end

    end

    
    methods
        function obj = Behavior2pPreprocessing(varargin)
            % obj = Behavior2pPreprocessing('ffm','C:\ffmpeg\bin\ffmpeg.exe','outdir','C:\rotated');
            %
            for i = 1:2:numel(varargin)
                k = varargin{i};
                v = varargin{i+1};
                if isprop(obj,k)
                    obj.(k) = v;
                end
            end
        end      
        function out = rotate_avis_ffmpeg(obj, src, deg, outdir, suf, ffm)
            % Rotate one or multiple AVI files using ffmpeg.
            %
            % src: folder, wildcard pattern, filename or cell array of filenames
            % deg: rotation angle in degrees
            % outdir (optional): output directory
            % suf (optional): suffix appended to output filename
            % ffm (optional): path to ffmpeg executable
            %
            % Returns:
            % out: string array of output file paths
        
            % Apply defaults from object properties
            if nargin < 4 || isempty(outdir)
                outdir = obj.outdir;
            end
            if nargin < 5 || isempty(suf)
                if ~isempty(obj.suf)
                    suf = obj.suf;
                else
                    suf = sprintf('_rot%gdeg',deg);
                end
            end
            if nargin < 6 || isempty(ffm)
                ffm = obj.ffm;
            end
        
            % Build list of input files from src
            if ischar(src) || isstring(src)
                src = char(src);
                if isfolder(src)
                    L = dir(fullfile(src,'*.avi'));
                    f = fullfile({L.folder},{L.name});
                elseif contains(src,'*')
                    L = dir(src);
                    f = fullfile({L.folder},{L.name});
                else
                    f = {src};
                end
            elseif iscell(src)
                f = src;
            else
                error('src must be folder, wildcard pattern, file, or cellstr');
            end
        
            % Output file path function
            if isempty(outdir)
                % Save next to original file
                g = @(p) fullfile(fileparts(p), ...
                    [erase(extname(p),'.avi'), suf, '.avi']);
            else
                % Save to outdir (create if needed)
                if ~exist(outdir,'dir')
                    mkdir(outdir);
                end
                g = @(p) fullfile(outdir, ...
                    [erase(extname(p),'.avi'), suf, '.avi']);
            end
        
            % Process each file
            nfiles = numel(f);
            out = cell(nfiles,1);
            parfor i = 1:nfiles
                inFile = f{i};
                outFile = g(inFile);
        
                % ffmpeg command with bilinear interpolation, MJPEG encoding
                cmd = sprintf(['%s -y -i "%s" ', ...
                    '-vf "rotate=%g*PI/180:bilinear=1" ', ...
                    '-c:v mjpeg -q:v %d -pix_fmt %s -an "%s"'], ...
                    ffm, inFile, deg, obj.q, obj.pix, outFile);
        
                % --- Verbose console output ---
                fprintf('[%d/%d] Rotating "%s" → "%s"\n',i,nfiles,inFile,outFile);
                fprintf('Running: %s\n',cmd);
        
                % Execute ffmpeg
                [status,cmdout] = system(cmd);
        
                % Print ffmpeg output if something goes wrong
                if status ~= 0
                    warning('ffmpeg failed for file "%s":\n%s',inFile,cmdout);
                end
        
                out{i} = string(outFile);
            end
        end
    end
end

% helpers
function e = extname(p)
% Return name+extension of file
[~,n,e2] = fileparts(p);
e = [n e2];
end

function s = quote(p)
%QU0TE Wrap a path or argument string in double quotes if not already quoted.
p = string(p);
if startsWith(p, '"') && endsWith(p, '"')
    s = char(p);
else
    s = char('"' + p + '"');
end

s = string(sprintf('"%s"', p));
end