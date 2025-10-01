classdef BinarizeLargeObjectOp < CropOp
    properties
        sens = 0          % adaptthresh sensitivity (0–1)
        win  = 101           % odd neighborhood size for adaptthresh
        sigma = 1.0          % Gaussian blur (px); 0 to disable
        min_area_frac = 0.002% remove blobs < this fraction of frame
        close_r = 3          % closing radius (px) to seal gaps; 0 to skip
        keep_largest = true  % keep only largest component
        profile = 'Grayscale AVI' % output profile
        out_suffix = '_mask'  % filename suffix for mask video
    end
    methods
        function o = BinarizeLargeObjectOp(varargin)
            for i = 1:2:numel(varargin)
                o.(varargin{i}) = varargin{i+1};
            end
        end

        function run(o, d, varargin)
            S = [dir(fullfile(d,'*.avi')); dir(fullfile(d,'*.mp4'))];
            for k = 1:numel(S)
                in = fullfile(S(k).folder, S(k).name);
                v = VideoReader(in);

                [~,n,~] = fileparts(in);
                out = fullfile(d, [n o.out_suffix '.avi']);
                w = VideoWriter(out, o.profile);
                w.FrameRate = v.FrameRate;
                open(w);
                c = onCleanup(@() try_close(w)); %#ok<NASGU>

                % prepare once from first frame
                f = readFrame(v);
                if size(f,3)==3, f = rgb2gray(f); end
                I = im2single(f);
                if o.sigma>0, I = imgaussfilt(I, o.sigma); end
                wn = max(3, o.win + mod(o.win+1,2)); % ensure odd
                T = adaptthresh(I, o.sens, 'NeighborhoodSize',[wn wn], 'ForegroundPolarity','bright');
                BW = imbinarize(I,T);
                Amins = ceil(numel(BW) * o.min_area_frac);

                BW = post(BW, o.close_r, o.keep_largest, Amins);
                writeVideo(w, uint8(BW)*255);

                while hasFrame(v)
                    f = readFrame(v);
                    if size(f,3)==3, f = rgb2gray(f); end
                    I = im2single(f);
                    if o.sigma>0, I = imgaussfilt(I, o.sigma); end
                    T = adaptthresh(I, o.sens, 'NeighborhoodSize',[wn wn], 'ForegroundPolarity','bright');
                    BW = imbinarize(I,T);
                    BW = post(BW, o.close_r, o.keep_largest, Amins);
                    writeVideo(w, uint8(BW)*255);
                end
            end
        end
    end
end

function BW = post(BW, close_r, keep_largest, Amin)
BW = bwareaopen(BW, Amin);                   % drop tiny specks
if keep_largest, BW = bwareafilt(BW,1); end  % keep big object
if close_r>0, BW = imclose(BW, strel('disk',close_r)); end
BW = imfill(BW,'holes');
end

function try_close(w)
try, close(w); end %#ok<TRYNC>
end
