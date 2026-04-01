classdef CropProject < handle
    % CropProject - manages named crops across large video shards
    %   - stores named ROI definitions
    %   - crops shards into per-name subfolders
    %   - dispatches on-demand operations (Strategy pattern)
    
    properties
        root              % project root folder
        crops             % Map: name -> struct('Xs','Ys','dir')
        ops               % Map: op name -> op object
        ang = 0           % rotation angle (deg CCW)
    end
    
    methods
        function o = CropProject(root)
            if nargin<1, root = pwd; end
            o.root = root;
            o.crops = containers.Map('KeyType','char','ValueType','any');
            o.ops   = containers.Map('KeyType','char','ValueType','any');
        end
        
        function set_rois(o, v_or_path, ang)
            % define crops by interactively drawing + naming
            if nargin<3 || isempty(ang), ang = 0; end
            o.ang = ang;
            if ischar(v_or_path) || isstring(v_or_path)
                v = VideoReader(v_or_path);
            else
                v = v_or_path;
            end
            m = defineLabeledBoxes(v,ang);
            ks = m.keys;
            for i=1:numel(ks)
                k = ks{i};
                d = fullfile(o.root, k);            % folder = name
                if ~exist(d,'dir'), mkdir(d); end
                s = m(k); s.dir = d;                % attach folder
                o.crops(k) = s;
            end
        end
        
        function crop_shards(o, shards)
            % crop all shards into per-crop folders
            if ischar(shards) || isstring(shards), shards = cellstr(shards); end
            for k=1:numel(shards)
                disp(['Cropping shard ', num2str(k), ' of ', num2str(numel(shards)), ': ', shards{k}]);
                o.crop_one_shard(shards{k})
            end
        end
        
        function crop_one_shard(o, shard_path)
            [~,n,~] = fileparts(shard_path);
            v = VideoReader(shard_path);
            ks = o.crops.keys;
            W = containers.Map('KeyType','char','ValueType','any');
            t = datestr(now,'yyyymmdd_HHMMSS');

            % choose container + profile
            ext = '.avi'; prof = 'Grayscale AVI';

            % create writers
            for i=1:numel(ks)
                k = ks{i}; d = o.crops(k).dir;
                %out = fullfile(d, sprintf('%s_%s_%s%s', n, k, t, ext));
                out = fullfile(d,  sprintf('%s%s', n(end-3:end), ext)); % idiosyncratic, brittle. reason: file path length
                w = VideoWriter(out, prof);
                w.FrameRate = v.FrameRate;
                open(w);
                W(k) = w;
            end

            % ensure closure even on error
            c = onCleanup(@() close_and_release(W));

            try
                h = waitbar(0, 'Processing frames...'); % Initialize progress bar
                frameCount = 0;
                while hasFrame(v)
                    f = readFrame(v);
                    frameCount = frameCount + 1;
                    if size(f,3)==3, f = rgb2gray(f); end
                    f = rot_frame(f, o.ang);
                    if ~isa(f,'uint8'), f = im2uint8(f); end

                    [H,Wid] = size(f);
                    for i=1:numel(ks)
                        k = ks{i}; bb = o.crops(k);
                        x1 = max(1, floor(bb.Xs(1))); x2 = min(Wid, ceil(bb.Xs(2))-1);
                        y1 = max(1, floor(bb.Ys(1))); y2 = min(H,   ceil(bb.Ys(2))-1);
                        if x2>x1 && y2>y1
                            writeVideo(W(k), f(y1:y2, x1:x2));
                        end
                    end
                    waitbar(frameCount / v.NumFrames, h); % Update progress bar
                end
                close(h); % Close progress bar
            catch ME
                rethrow(ME)  % onCleanup will still run
            end
        end
        
        function d = crop_dir(o, name)
            % return folder path of a named crop
            d = o.crops(name).dir;
        end
        
        function names = crop_names(o)
            % return all crop names
            names = o.crops.keys;
        end
        
        function register_op(o, name, op_obj)
            % store operation object under a key
            arguments
                o 
                name (1,:) {mustBeText}
                op_obj (1,1) CropOp
            end

            o.ops(name) = op_obj;
        end
        
        function run_op(o, op_name, crop_name, varargin)
            % run one op on one named crop

            % this formulation was necessary to avoid conflict with
            % shadowed names
            op = o.ops(op_name);
            d  = o.crops(crop_name).dir;
            op.run(d,varargin{:})

            % o.ops(op_name).run(o.crop_dir(crop_name), varargin{:});
        end
        
        function run_op_all(o, op_name, varargin)
            % run op on all crops
            ks = o.crops.keys;
            for i=1:numel(ks), o.run_op(op_name, ks{i}, varargin{:}); end
        end
        
        function L = list_crop_files(~, dirpath, ext)
            % list files inside a crop folder
            if nargin<3, ext='*.avi'; end
            s = dir(fullfile(dirpath,ext));
            L = fullfile({s.folder},{s.name});
        end
    end
end

%% Helpers

function close_and_release(W)
    if isempty(W), return; end
    ks = W.keys;
    for i=1:numel(ks)
        w = W(ks{i});
        try, close(w); end %#ok<TRYNC>
    end
    for i=1:numel(ks)      % drop references to unlock files on Windows
        remove(W, ks{i});
    end
    clear W
    pause(0.05)            % tiny yield to let OS release handles
end

function g = rot_frame(f, ang)
    a = mod(ang,360);
    r = round(a/90);
    if abs(a-90*r)<1e-9, g = rot90(f, r);
    else, g = imrotate(f, ang, 'nearest', 'crop');
    end
end