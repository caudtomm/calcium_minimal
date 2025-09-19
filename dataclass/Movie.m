classdef Movie
    % A Movie Object
    properties
        description string = ""
        detailed_description string = ""

        stack
        path
        timeavg double
        badperiods double = []
        log cell = {}

        % initially copied from scanimage_meta for quicker access
        h double = 0    % frame height
        w double = 0    % frame width
        nfr double = 1  % number of frames
        scanimage_meta cell
        metadata struct

    end
    properties (SetAccess = protected)
        fs double       % framerate
        frameperiod double = 0
        imagingFovUm double = [0 0; 0 0; 0 0; 0 0]  % frame size [um]
    end

    methods
        function obj = Movie(src) % Constructor
            % arg 'src' is filename or numerical matrix
            arguments
                src
            end

            if nargin>0 && ~isempty(src)
                [obj.stack, obj.scanimage_meta, obj.path, meta] = readSource(src);
                % populate properties from meta when present (e.g., AVI case)
                if ~isempty(meta)
                    if isfield(meta,'fs') && ~isempty(meta.fs)
                        obj.fs = meta.fs;
                        obj.frameperiod = 1/meta.fs;
                    end
                    if isfield(meta,'frameperiod') && ~isempty(meta.frameperiod)
                        obj.frameperiod = meta.frameperiod;
                        if isempty(obj.fs) || obj.fs==0
                            obj.fs = 1/meta.frameperiod;
                        end
                    end

                    obj.metadata = meta;
                end
            end

            % copy properties from 'scanimage_meta'
            try
                obj.h = obj.scanimage_meta{1}.SI_hRoiManager_linesPerFrame;
                obj.w = obj.scanimage_meta{1}.SI_hRoiManager_pixelsPerLine;
                obj.nfr = numel(obj.scanimage_meta);
                obj.fs = obj.scanimage_meta{1}.SI_hRoiManager_scanFrameRate;
                obj.frameperiod = obj.scanimage_meta{1}.SI_hRoiManager_scanFramePeriod;
                obj.imagingFovUm = obj.scanimage_meta{1}.SI_hRoiManager_imagingFovUm;
            catch
                % fallback if scanimage_meta is empty or missing fields
                [obj.h,obj.w,obj.nfr] = size(obj.stack);
            end
        end

        % getters
        function value = get.timeavg(obj)
            if (isa(obj.stack,'StackH5') || isa(obj.stack,'StackH5Multi'))
                if isfield(obj.metadata,'timeavg_image') && ~isempty(obj.metadata.timeavg_image)
                    value = obj.metadata.timeavg_image;    % O(1): cached
                    return
                end
                value = mean(obj.stack,3,'omitmissing');    % fallback: stream compute
                return
            end
            value = mean(obj.stack,3,'omitmissing');        % in-RAM double
        end

        function value = get.h(obj)
            if isfield(obj.metadata,'h') && ~isempty(obj.metadata.h)
                value = obj.metadata.h;                     % O(1): cached
                return
            end
            value = size(obj.stack,1);                      % cheap even for proxies
        end

        function value = get.w(obj)
            if isfield(obj.metadata,'w') && ~isempty(obj.metadata.w)
                value = obj.metadata.w;                     % O(1): cached
                return
            end
            value = size(obj.stack,2);
        end

        function value = get.nfr(obj)
            if isfield(obj.metadata,'nfr') && ~isempty(obj.metadata.nfr)
                value = obj.metadata.nfr;                   % O(1): cached
                return
            end
            value = size(obj.stack,3);
        end

        function obj = setFramerateHz(obj, value)
            obj.fs = value;
            obj.frameperiod = 1/value;
        end

        function obj = setFramePeriodSec(obj, value)
            obj = obj.setFramerateHz(1/value);
        end

        function obj = fliplr(obj)
            obj.stack = fliplr(obj.stack);
            obj = obj.checkAgain;
        end

        function obj = flipud(obj)
            obj.stack = flipud(obj.stack);
            obj = obj.checkAgain;
        end

        function obj = checkAgain(obj) %% EMPTY: TESTING IF NO LONGER NECESSARY (REPLACED BY MORE EFFICIENT GETTERS)
            % %initialize
            % movie = obj.stack;
            % 
            % % re-compute
            % [height,width,framecount] = size(movie);
            % zproj = mean(movie,3,'omitmissing');
            % 
            % % store
            % obj.h = height;
            % obj.w = width;
            % obj.nfr = framecount;
            % obj.timeavg = zproj;
        end

        function v = getmaxval(obj)
            if isfield(obj.metadata,'mean_all') && isfield(obj.metadata,'std_all') ...
                    && ~isempty(obj.metadata.mean_all) && ~isempty(obj.metadata.std_all)
                v = obj.metadata.mean_all + 4*obj.metadata.std_all; return
            end
            v = mean(obj.stack,"all","omitmissing") + 4*std(obj.stack,[],"all","omitmissing");
        end

        
        % Method to log operations and their hash address
        function obj = update_log(obj,charv)
            % append a new voice to log
            obj.log{end+1} = charv;
            % make log vertical for easier visualization (only needs to happen once)
            if numel(obj.log)==2; obj.log = transpose(obj.log); end
        end

        % function to play the movie
        function play(obj,framerate_Hz)
            arguments
                obj
                framerate_Hz double = obj.fs
            end
            if isempty(framerate_Hz) || ~isfinite(framerate_Hz) || framerate_Hz<=0
                framerate_Hz = 16;
            end
            dt = 1/framerate_Hz;

            % Fixed display spans (avoid per-frame min/max or autoscale)
            % If your data are uint8, this is perfect. Otherwise keep getmaxval.
            if isa(obj.stack,'StackH5') || isa(obj.stack,'StackH5Multi')
                cspan = [0 255];   % disk-backed path stores uint8
                xspan = [0 255];
            else
                cspan = [0 obj.getmaxval];
                xspan = [0 obj.getmaxval];
            end

            % Figure & layout
            hfig = figure('Name','Movie','NumberTitle','off','Color','w');
            setappdata(hfig,'stop',false);
            hfig.CloseRequestFcn = @(src,evt) setappdata(src,'stop',true);  % mark stop; we'll delete later
            hfig.Position = [50 50 1500 800];
            tl = tiledlayout(hfig,1,2,'Padding','compact','TileSpacing','compact');
            colormap(hfig,"gray")

            % First frame
            f = 1;
            frame = obj.stack(:,:,f);                 % should be uint8 for speed

            % Left: image (use image + fixed CLim)
            ax1 = nexttile(tl,1);
            him = image(ax1, frame, 'CDataMapping','scaled');
            ax1.CLim = cspan;
            axis(ax1,'image','off')
            % Single text label instead of title (faster)
            htxt = text(ax1, 10, 20, sprintf('frame: %d',f), 'Color','w', ...
                'FontWeight','bold', 'FontName','Helvetica', 'FontSize',12);

            % Right: histogram (outline via stairs OR area; here we use area, throttled)
            ax2 = nexttile(tl,2);
            nbins = 256;
            edges = linspace(xspan(1), xspan(2), nbins+1);
            centers = (edges(1:end-1)+edges(2:end))/2;
            % imhist is fast for uint8 (IPT), else use histcounts
            if isa(frame,'uint8')
                counts = imhist(frame, nbins);
                % imhist returns 256×1 for uint8; centers must match 1:256 mapping
                centers = linspace(0,255,nbins);
            else
                counts = histcounts(frame(:), edges);
            end
            harea = area(ax2, centers, counts);
            harea.LineWidth = 0.5; harea.EdgeAlpha = 0.5; harea.FaceAlpha = 0.25;
            xlim(ax2, xspan); ax2.YLimMode = 'auto'; ax2.Box = 'on';

            hold(ax2,'on')
            l1 = line(ax2,[cspan(1) cspan(1)],[0 1],'Color','r','LineStyle','--','LineWidth',2);
            l2 = line(ax2,[cspan(2) cspan(2)],[0 1],'Color','r','LineStyle','--','LineWidth',2);
            hold(ax2,'off')

            % Playback loop with frame skipping + throttled histogram updates
            HIST_EVERY = 3;                 % update histogram every N frames
            next_hist  = 1;
            t0 = tic;
            frame_start_time = 0;           % target time for frame 1 (sec)

            while isvalid(hfig)
                for f = 1:obj.nfr
                    if ~ishghandle(hfig) || getappdata(hfig,'stop'); safe_close(hfig); return; end
                    if ~isvalid(hfig); return; end

                    % Catch up if we're behind: skip frames
                    now = toc(t0);
                    target = frame_start_time + (f-1)*dt;
                    if now > target + dt
                        % Skip ahead by however many frames we're behind
                        behind = floor((now - target)/dt);
                        f = min(obj.nfr, f + behind);
                    end

                    % Grab frame (fast path returns 2D uint8)
                    frame = obj.stack(:,:,f);

                    % Update image + label
                    set(him,'CData',frame);
                    set(htxt,'String',sprintf('frame: %d',f));

                    % Throttled histogram update
                    if mod(f, HIST_EVERY)==next_hist
                        if isa(frame,'uint8')
                            counts = imhist(frame, nbins);
                            centers = linspace(0,255,nbins);
                        else
                            counts = histcounts(frame(:), edges);
                        end
                        set(harea,'XData',centers,'YData',counts);
                        ymax = max(counts); if ~isfinite(ymax) || ymax<=0, ymax = 1; end
                        set(l1,'YData',[0 ymax*2/3]); set(l2,'YData',[0 ymax*2/3]);
                    end

                    drawnow limitrate

                    % Sleep the remainder to hit target dt (if any)
                    now = toc(t0);
                    to_sleep = target + dt - now;
                    if to_sleep > 0
                        pause(to_sleep);
                    end
                end
                % loop continuously
                frame_start_time = toc(t0);
            end
        end

        function playFast(obj, framerate_Hz, useGPU, prefetchN)
            arguments
                obj
                framerate_Hz double = obj.fs
                useGPU logical = true
                prefetchN double {mustBeInteger, mustBePositive} = 32   % frames buffered ahead
            end

            if isempty(framerate_Hz) || ~isfinite(framerate_Hz) || framerate_Hz<=0
                framerate_Hz = 16;
            end
            dt = 1/framerate_Hz;

            % Ensure we use the GPU for figure rendering if available
            try, opengl hardware; catch, end

            % Choose display range & type (uint8 path is fastest)
            isDiskBacked = isa(obj.stack,'StackH5') || isa(obj.stack,'StackH5Multi');
            if isDiskBacked
                cspan = [0 255]; xspan = [0 255]; toUint8 = @(x) x;   % shards are uint8
            else
                cspan = [0 obj.getmaxval]; xspan = cspan;
                toUint8 = @(x) im2uint8(mat2gray(x, cspan));           % convert once per frame
            end

            % Figure (lean)
            hfig = figure('Name','Movie','NumberTitle','off','Color','w', ...
                        'MenuBar','none','ToolBar','none','Renderer','opengl');
            setappdata(hfig,'stop',false);
            hfig.CloseRequestFcn = @(src,evt) setappdata(src,'stop',true);
            hfig.Position = [50 50 1500 800];
            set(hfig,'GraphicsSmoothing','off');

            tl = tiledlayout(hfig,1,2,'Padding','compact','TileSpacing','compact');
            colormap(hfig, gray(256));  % fixed LUT for indexed display

            % First frame (synchronously)
            f = 1;
            frame = obj.stack(:,:,f);                   % expect uint8 for disk-backed path
            frameU8 = toUint8(frame);

            % LEFT: show as indexed image (fastest path)
            ax1 = nexttile(tl,1);
            him = image(ax1, frameU8, 'CDataMapping','direct');  % uint8 indices -> colormap
            axis(ax1,'image','off')
            ax1.CLim = [0 255];                         % fixed
            htxt = text(ax1,10,20,sprintf('frame: %d',f),'Color','w','FontWeight','bold');

            % RIGHT: histogram (GPU-optional, throttled)
            ax2 = nexttile(tl,2);
            nbins = 256;
            centers = uint8(0:255);
            if useGPU && canUseGPU(); 
                counts = gather(histcounts(gpuArray(frameU8), nbins, 'BinLimits',[0,255]));
            else
                counts = histcounts(frameU8(:), nbins, 'BinLimits',[0,255]);
            end
            harea = area(ax2, double(centers), counts);
            harea.LineWidth = 0.5; harea.EdgeAlpha = 0.5; harea.FaceAlpha = 0.25;
            xlim(ax2, xspan); ax2.YLimMode = 'auto'; ax2.Box = 'on';
            hold(ax2,'on')
            l1 = line(ax2,[cspan(1) cspan(1)],[0 1],'Color','r','LineStyle','--','LineWidth',2);
            l2 = line(ax2,[cspan(2) cspan(2)],[0 1],'Color','r','LineStyle','--','LineWidth',2);
            hold(ax2,'off')

            % --- Prefetch setup (background worker) ---
            import parallel.pool.*
            dq = PollableDataQueue;      % main thread polls frames
            doneFlag = parallel.pool.DataQueue; % worker reports finish/errors
            setappdata(hfig,'dq_done',false);
            afterEach(doneFlag, @(msg) setappdata(hfig,'dq_done',true));

            F = obj.nfr;
            % Start worker to prefetch in batches
            if F>1
                pool = gcp('nocreate');
                if isempty(pool), parpool('threads'); end  % light-weight pool
                % Pass a lightweight handle to the reader via function handle
                fcn = @() prefetch_worker(obj, dq, doneFlag, prefetchN);
                parfeval(@() fcn(), 0);   % fire and forget
            end

            % Playback loop with frame skipping + throttled histogram
            HIST_EVERY = 4;  next_hist = 1;
            t0 = tic; frame_start_time = 0;
            f = 1;

            while isvalid(hfig)
                if getappdata(hfig,'stop'), safe_close(hfig); return; end

                % Try to poll a prefetched frame; fallback to synchronous read
                data = poll(dq, 0);   % non-blocking
                if ~isempty(data)
                    f = data.frameIdx;
                    frameU8 = data.frameU8;
                else
                    % synchronous (should be rare if prefetch keeps up)
                    frameU8 = toUint8(obj.stack(:,:,f));
                end

                % Catch up if behind: skip forward
                now = toc(t0); target = frame_start_time + (f-1)*dt;
                if now > target + dt
                    behind = floor((now - target)/dt);
                    f = min(F, f + behind);
                end

                % Update visuals
                set(him,'CData',frameU8);
                set(htxt,'String',sprintf('frame: %d',f));

                if mod(f, HIST_EVERY)==next_hist
                    if useGPU && canUseGPU()
                        counts = gather(histcounts(gpuArray(frameU8), nbins, 'BinLimits',[0,255]));
                    else
                        counts = histcounts(frameU8(:), nbins, 'BinLimits',[0,255]);
                    end
                    set(harea,'XData',double(centers),'YData',counts);
                    ymax = max(counts); if ~isfinite(ymax) || ymax<=0, ymax = 1; end
                    set(l1,'YData',[0 ymax*2/3]); set(l2,'YData',[0 ymax*2/3]);
                end

                drawnow limitrate   % allow callbacks (close button)
                if ~ishghandle(hfig) || getappdata(hfig,'stop'), safe_close(hfig); return; end

                % Sleep remainder to hit FPS
                now = toc(t0);
                to_sleep = target + dt - now;
                if to_sleep > 0, pause(to_sleep); end

                % Next frame (wrap)
                f = f + 1; if f>F, f = 1; frame_start_time = toc(t0); end
            end
        end

        function FileOut = save(obj, newpath, type, newfname, auto_overwrite)
            arguments
                obj Movie
                newpath char = ''
                type char = 'mat'
                newfname char = ''
                auto_overwrite logical = false
            end
            
            outpath = pwd;
            if ~isempty(newpath)
                outpath = newpath;
            % elseif ~isempty(obj.path) && ~isempty(obj.path.orig_fpath)
            %     outpath = obj.path.orig_fpath;
            end

            outfname = 'movie';
            if ~isempty(newfname)
                outfname = newfname;
            elseif ~isempty(obj.path) && ~isempty(obj.path.fname)
                outfname = obj.path.fname;
            end

            % check if the output path exists. if not, create it.
            if ~exist(outpath,"dir"); mkdir(outpath); end

            switch type
                case 'mat'
                    FileOut = fullfiletol(outpath,[outfname,'.mat']);
                    b = true;
                    if ~auto_overwrite
                        b = prompt_overwrite(FileOut);
                    end
        
                    obj.path = getFileNameSpecs(FileOut);
        
                    movie = obj;
                    if b
                        s.movie = movie;
                        robust_io('save',FileOut,s);
                    end
                    
                case 'tif'
                    FileOut = fullfiletol(outpath,[outfname,'.tif']);
                    totpathlen = length(fullfiletol(pwd,FileOut));
                    if totpathlen>259
                        if ~isempty(getFileNameSpecs(outfname).trial_num)
                            outfname = ['tr_',num2str(getFileNameSpecs(outfname).trial_num)];
                        else
                            outfname = 'filepath_too_long';
                            % # TODO: This should call a little function to 
                            % pick a unique filename.
                            % Also, this whole if-block should be extruded
                            % from the switch block and applied in common
                            % between classes Movie and Snippet (adapt use
                            % suitably for each). Or maybe it should be its
                            % own function and called for each 'type'case :
                            % more flexibility for extensions of different
                            % length?
                        end
                        FileOut = fullfiletol(outpath,[outfname,'.tif']);
                        totpathlen2 = length(fullfiletol(pwd,FileOut));
                        fprintf('\nOutput file path is too long! Lenght: %s > 259\nNew filename: %s\nNew length: %s\n\n', ...
                            num2str(totpathlen), outfname, num2str(totpathlen2));
                    end
                    b = true;
                    if ~auto_overwrite
                        b = prompt_overwrite(FileOut);
                    end
        
                    if b
                        delete(FileOut)
                        if isa(obj.stack,'StackH5')
                            % stream-write HDF5-backed stack to TIFF without loading all frames
                            t = Tiff(FileOut,'w8');
                            tag.ImageLength = obj.h; tag.ImageWidth = obj.w;
                            tag.Compression = Tiff.Compression.None;
                            tag.Photometric = Tiff.Photometric.MinIsBlack;
                            tag.BitsPerSample = 16; tag.SamplesPerPixel = 1;
                            tag.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
                            for frame = 1:obj.nfr
                                setTag(t,tag);
                                write(t, uint16(obj.stack(:,:,frame)));
                                if frame<obj.nfr, writeDirectory(t); end
                            end
                            close(t)
                        else
                            movie = obj.stack;
                            saveastiff(uint16(movie), FileOut)
                        end
                    end

                case 'avi'
                    FileOut = fullfiletol(outpath,[outfname,'.avi']);
                    b = true;
                    if ~auto_overwrite
                        b = prompt_overwrite(FileOut);
                    end

                    if b
                        v = VideoWriter(FileOut);
                        v.FrameRate = obj.fs;
                        open(v);

                        for frame = 1:obj.nfr
                            normframe = obj.stack(:, :, frame)/obj.getmaxval;
                            normframe(normframe>1) = 1;
                            normframe(normframe<0) = 0;
                            writeVideo(v,normframe);
                        end
                        close(v)
                    end


                otherwise
                    error(['File type ''',type,''' not supported for saving.'])
            end
        end
    end
end


function [stack,scanimage_meta,path,meta] = readSource(src)
    % initialize
    scanimage_meta = {};
    [stack, path] = deal([]);
    meta = struct(); % optional non-ScanImage metadata (e.g., AVI info)

    % read
    if isnumeric(src)
        src = double(src);
    end
    switch class(src)
        case {'string','char'}  % src :
            src = char(src);
            disp(['Loading: ',src])
            if endsWith(src,'.mat')
                movie = robust_io('load',src,'movie').movie;
                stack = movie.stack;
                scanimage_meta = movie.scanimage_meta;
                path = movie.path;
                try
                    if isprop(movie,'fs') && ~isempty(movie.fs)
                        meta.fs = movie.fs;
                        meta.frameperiod = 1/movie.fs;
                    end
                    if isprop(movie,'h') && ~isempty(movie.h), meta.h = movie.h; end
                    if isprop(movie,'w') && ~isempty(movie.w), meta.w = movie.w; end
                    if isprop(movie,'nfr') && ~isempty(movie.nfr), meta.nfr = movie.nfr; end
                catch
                end
            elseif endsWith(src,'.tif') || endsWith(src,'.tiff')
                [stack, scanimage_meta] = loadTiffStack(src);
                stack = double(stack);
                try
                    scanimage_meta = convertScanimageMeta(scanimage_meta);
                catch
                end
                path = getFileNameSpecs(src);
                        elseif endsWith(src,'.avi')
                [stack, meta] = ingestAvi(src, 10); % GB limit for double
                path = getFileNameSpecs(src);
            end
            if isempty(path.orig_fpath)
                path.orig_fpath = pwd;
            end
        case {'double'}         % src : numeric matrix defining the movie
            stack = src;
            path.fname = 'movie';
    end
end

function result = convertScanimageMeta(char_meta)
    result = char_meta;
    for i = 1:numel(char_meta)
        result{i} = parseToStruct(char_meta{i});
    end
end

function result = parseToStruct(inputStr)
    % Initialize an empty structure
    result = struct();
    
    % Split the input string into lines
    lines = strsplit(inputStr, '\n');
    
    % Process each line
    for i = 1:length(lines)
        line = strtrim(lines{i});
        if isempty(line)
            continue;
        end
        
        % Find the equals sign to split parameter names and values
        eqIndex = find(line == '=', 1);
        if isempty(eqIndex)
            continue;
        end
        
        % Extract the parameter name and value
        paramName = strtrim(line(1:eqIndex-1));
        paramValue = strtrim(line(eqIndex+1:end));
        
        % Check if the value is numeric, a string, or an array
        if ~isempty(paramValue) && all(ismember(paramValue(1), '0123456789+-.[],'))
            % Try to evaluate the expression to handle numbers and vectors
            try
                paramValue = eval(paramValue);
            catch
                % If eval fails, revert to treating it as a string
                paramValue = strrep(paramValue, '''', ''); % Remove any enclosing single quotes
            end
        else
            % It's a string, strip single quotes
            paramValue = strrep(paramValue, '''', '');
        end
        
        % Assign the value to the structure
        paramName = matlab.lang.makeValidName(paramName); % Ensure the field name is a valid MATLAB identifier
        result.(paramName) = paramValue;
    end
end

function [stk,meta] = ingestAvi(src,limGB)
v = VideoReader(src);
F = floor(v.Duration*v.FrameRate);
H = v.Height; W = v.Width; fs = v.FrameRate;
est = F*H*W*8;
if est > limGB*1e9
    [files,lens,S,sum1,sum2,N] = avi_shard_write(src,H,W,F);
    stk = StackH5Multi(files,lens,[H W F]);
    mu = sum1/N;
    sd = sqrt(max((sum2 - (sum1^2)/N)/max(N-1,1),0));
    meta.backend='hdf5_shards';
    meta.shards=files; meta.shard_len=lens;
    meta.timeavg_image = S/F; meta.mean_all=mu; meta.std_all=sd;
    meta.bitdepth='uint8'; meta.colorspace='grayscale';
else
    stk = zeros(H,W,F,'double');
    bp = v.BitsPerPixel/3; sc = 2^bp-1;
    k=1; while hasFrame(v)
        f = readFrame(v);
        g = rgb2gray(f/sc)*sc;
        stk(:,:,k) = double(g); k=k+1;
    end
    meta.backend='memory'; meta.bitdepth='double_from_uint8'; meta.colorspace='grayscale';
    meta.timeavg_image = mean(stk,3);
    meta.mean_all = mean(stk(:));
    meta.std_all  = std(stk(:),0);
end
meta.source='avi'; meta.h=H; meta.w=W; meta.nfr=F; meta.fs=fs; meta.frameperiod=1/fs;
end

function [files,lens,S,sum1,sum2,N] = avi_shard_write(src,H,W,F)
od = fullfile(tempdir, 'movie_h5_shards');
if ~exist(od,'dir'), mkdir(od); end
d = dir(fullfile(od,'shard_*.h5')); for i=1:numel(d), try, delete(fullfile(od,d(i).name)); end, end

p = gcp('nocreate'); n = 1; if ~isempty(p), n = max(1,p.NumWorkers); end
n = min(n, max(1, ceil(F/256)));
cuts = round(linspace(1,F+1,n+1));
rngs = [cuts(1:end-1).' (cuts(2:end)-1).'];

files = strings(size(rngs,1),1);
lens  = zeros(size(rngs,1),1);

S_cell = cell(size(rngs,1),1);
s1 = zeros(size(rngs,1),1);
s2 = zeros(size(rngs,1),1);
nn = zeros(size(rngs,1),1);

B = 128; % frames per write slab

parfor i = 1:size(rngs,1)
    a = rngs(i,1); b = rngs(i,2); L = b-a+1;
    v = VideoReader(src); v.CurrentTime = (a-1)/v.FrameRate;

    hf = fullfile(od,sprintf('shard_%03d.h5',i));
    h5create(hf,'/stack',[H W L],'Datatype','uint8','Chunksize',[H W min(L,128)],'Deflate',0);

    S_loc = zeros(H,W,'double');
    s1_loc = 0; s2_loc = 0; n_loc = 0;

    off = 1;
    while off <= L
        nwrite = min(B, L - off + 1);
        buf = zeros(H,W,nwrite,'uint8');
        for k = 1:nwrite
            f = readFrame(v);
            buf(:,:,k) = rgb2gray(f);
        end
        h5write(hf,'/stack',buf,[1 1 off],[H W nwrite]);

        bd = double(buf);
        S_loc = S_loc + sum(bd,3);
        s1_loc = s1_loc + sum(bd(:));
        s2_loc = s2_loc + sum(bd(:).^2);
        n_loc  = n_loc  + numel(bd);

        off = off + nwrite;
    end

    files(i) = string(hf);
    lens(i)  = L;
    S_cell{i} = S_loc;
    s1(i) = s1_loc; s2(i) = s2_loc; nn(i) = n_loc;
end

S = zeros(H,W,'double');
for i=1:numel(S_cell), S = S + S_cell{i}; end
sum1 = sum(s1); sum2 = sum(s2); N = sum(nn);
end

function prefetch_worker(obj, dq, doneFlag, prefetchN)
    % simple ring buffer prefetcher
    try
        F = obj.nfr; k = 1;
        while true
            % fill up to prefetchN frames ahead
            for j = 1:prefetchN
                fidx = k + j - 1; if fidx>F, fidx = fidx - F; end
                fr = obj.stack(:,:,fidx);     % should be uint8 (fast path)
                send(dq, struct('frameIdx', fidx, 'frameU8', fr));
            end
            k = k + prefetchN;
            if k>F, k = k - F; end
            pause(0.001); % yield a bit
        end
    catch
        send(doneFlag, true);
    end
end

function tf = canUseGPU()
    try, tf = parallel.gpu.GPUDevice.isAvailable; catch, tf = false; end
end

function safe_close(h)
    if ishghandle(h), try, delete(h); end, end
end