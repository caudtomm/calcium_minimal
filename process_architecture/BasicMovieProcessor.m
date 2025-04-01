classdef BasicMovieProcessor < MovieProcessing
% Implements small functions for processing movies.
% Processors available:
%   - 'subtract_baseline'
%   - 'remove_badperiods'
%   - 'add_badperiods_as_nans'
%   - 'clahe'
%   - 'movmean'
%   - 'downsamplet'
%   - 'gauss_blur2d'
%
    
    properties (SetAccess = private)
        processor char
    end
    properties
        data_processed
        operation
    end

    methods
        function obj = BasicMovieProcessor(processor,movie_in)
            arguments
                processor char
                movie_in = ''
            end
            obj = obj@MovieProcessing(movie_in);
            if ~isempty(processor)
                obj = obj.setProcessor(processor);
            end

            % initialization
            obj.data_processed = obj.data_raw;
        end

        function obj = setProcessor(obj,processor)
            arguments
                obj BasicMovieProcessor
                processor char
            end
            obj.processor = processor;
            obj.init.method = processor;
        end

        function obj = run(obj,varargin)
            disp(['Running: ', obj.processor])
            switch obj.processor

                case 'subtract_baseline'
                    [movie, obj.operation.baseline] = ...
                        subtract_movie_baseline(obj.data_raw);
                
                case 'remove_badperiods'
                    [movie, obj.operation.periods_removed] = ...
                        subtract_badperiods(obj.data_raw);
                
                case 'add_badperiods_as_nans'
                    [movie, obj.operation.periods_added] = ...
                        add_badperiods_as_nans(obj.data_raw);

                case 'replace_badperiods_with_nans'
                    [movie, obj.operation.periods_removed] = ...
                        subtract_badperiods(obj.data_raw);
                    [movie, obj.operation.periods_added] = ...
                        add_badperiods_as_nans(movie.checkAgain);
                
                case 'clahe'
                    disp('Histogram equalization (CLAHE)')
                    [movie, obj.operation] = clahe(obj.data_raw);

                case 'movmean'
                    % Default window size
                    try
                        win_size = obj.operation.win_size;
                    catch
                        win_size = 5; 
                    end

                    % Parse varargin for win_size
                    if ~isempty(varargin)
                        for i = 1:length(varargin)-1
                            if ischar(varargin{i}) && strcmp(varargin{i}, 'win_size')
                                win_size = varargin{i+1};
                            end
                        end
                    end

                    [movie, obj.operation.win_size] = ...
                        movie_movmean(obj.data_raw, win_size);

                case 'downsamplet'
                    % Default downsampling factor
                    try
                        factor = obj.operation.factor;
                    catch
                        factor = 2; 
                    end

                    % Parse varargin for factor
                    if ~isempty(varargin)
                        for i = 1:length(varargin)-1
                            if ischar(varargin{i}) && strcmp(varargin{i}, 'factor')
                                factor = varargin{i+1};
                            end
                        end
                    end
                    [movie, obj.operation.factor] = ...
                        movie_downsample_time(obj.data_raw, factor);

                case 'gauss_blur2d'
                    % Default kernel sigma
                    try
                        sigma = obj.operation.sigma;
                    catch
                        sigma = 5; 
                    end

                    % Parse varargin for sigma
                    if ~isempty(varargin)
                        for i = 1:length(varargin)-1
                            if ischar(varargin{i}) && strcmp(varargin{i}, 'sigma')
                                sigma = varargin{i+1};
                            end
                        end
                    end
                    [movie, obj.operation.sigma] = ...
                        gauss_blur2d(obj.data_raw, sigma);

                otherwise
                    error('Basic movie processor function unknown.')
            end

            obj.data_processed = movie;
            obj = obj.tailsequence(obj.processor);
        end
    end
end

function [movie_subtracted, baseline] = subtract_movie_baseline(data_raw)
arguments
    data_raw Movie
end
movie = data_raw.stack;
framerate = data_raw.fs;
movie_subtracted = data_raw;
% baseline defined as the 3% quantile of all pixel values per frame (50 sec rolling avg)
baseline = movmean(quantile(movie,0.03,[1,2]),floor(50*framerate),'omitnan'); 

baselineMAT = repmat(baseline,[size(movie,1),size(movie,2)]);
movie_subtracted.stack = movie - baselineMAT;

baseline = squeeze(baseline);
end

function [movie_without_bp, adapted_badperiods] = subtract_badperiods(data_raw)
% Subtract bad periods
arguments
    data_raw Movie
end

% initialize
movie_without_bp = data_raw;
movie = data_raw.stack;
badperiods = data_raw.badperiods;
adapted_badperiods = [];
if isempty(badperiods); return; end
numberframes = data_raw.nfr;

% adapt bad periods to a sub-trial movie
for i_bp = 1:size(badperiods,1)
    thisbp = badperiods(i_bp,:);
    if thisbp(2)>numberframes; continue; end
    if thisbp(3)>numberframes; thisbp(3)=numberframes; end
    adapted_badperiods = [adapted_badperiods; thisbp];
end

for i_bp = 1:size(badperiods,1)
    thisbp = adapted_badperiods(size(adapted_badperiods,1)-i_bp+1,:);
    movie(:,:,thisbp(2):thisbp(3)) = [];
end
movie_without_bp.stack = movie;

fprintf('Removed following periods:\n');
fprintf('%g\t%g\t%g\n', adapted_badperiods.');
fprintf('\n');
disp(['New frame count: ',num2str(data_raw.nfr), ' --> ', ...
    num2str(size(movie_without_bp.stack,3)),' fr']);
end

function [movie_result, badperiods] = add_badperiods_as_nans(data_raw)
% add nans in the place of bad periods
arguments
    data_raw Movie
end

% initialize
movie_result = data_raw;
movie = data_raw.stack;
badperiods = data_raw.badperiods;
if isempty(badperiods); return; end
h = data_raw.h;
w = data_raw.w;
numbadframes = sum(diff(badperiods(:,2:3),[],2)+1);
nfr = data_raw.nfr+numbadframes;

if isempty(badperiods);return;end

tmp = nan(h,w,nfr);
bptab = [0 0 0; badperiods; 0 nfr+1 0];
for i_bp = 1:size(badperiods,1)+1
    idx = sum(bptab(1:i_bp+1,2)) - sum(bptab(1:i_bp,3));
    dur = bptab(i_bp+1,2)-bptab(i_bp,3)-1;
    tmp(:,:,bptab(i_bp,3)+1:bptab(i_bp+1,2)-1) = ...
        movie(:,:,idx-dur-i_bp+1:idx-i_bp);
end
movie_result.stack = tmp;

fprintf('Added following periods:\n');
fprintf('%g\t%g\t%g\n', badperiods.');
fprintf('\n');
disp(['New frame count: ',num2str(data_raw.nfr), ' --> ', ...
    num2str(size(movie_result.stack,3)),' fr']);
end

function [movie_result, settings] = clahe(data_raw)
arguments
    data_raw Movie
end

% initialize
movie_result = data_raw;
movie = data_raw.stack;
framecount = data_raw.nfr;

% hardcoded settings
settings.clipLimit = 0.01;
settings.Distribution = 'uniform';
settings.Range = 'original';
settings.NumTiles = [2 2];

% run
for i_f = 1:framecount
    movie(:,:,i_f) = adapthisteq(mat2gray(movie(:,:,i_f)), ...
        'clipLimit',settings.clipLimit, ...
        'Distribution',settings.Distribution);
end

% store to results
movie_result.stack = movie;

end

function [movie_result, win_size] = movie_movmean(data_raw, win_size)
arguments
    data_raw Movie
    win_size double
end

fprintf('Averaging over window size: %s frames\n',num2str(win_size))
    
% initialize
movie_result = data_raw;
movie = data_raw.stack;

movie = permute(movmean(permute(movie,[3,1,2]),win_size),[2,3,1]);

% store to results
movie_result.stack = movie;

end

function [movie_result, factor] = movie_downsample_time(data_raw, factor)
arguments
    data_raw Movie
    factor double
end

newfs = data_raw.fs/factor;

fprintf('Downsampling in time by a factor of %s - FR %s -> %sHz\n',...
    num2str(factor),num2str(data_raw.fs),num2str(newfs))
    
% initialize
movie_result = data_raw;
movie = data_raw.stack;

% execute
idx = mod([1:data_raw.nfr],factor) ~= 0;
movie(:,:,idx) = [];

% store to results
movie_result.stack = movie;
movie_result = movie_result.setFramerateHz(newfs);

end

function [movie_result, sigma] = gauss_blur2d(data_raw, sigma)
arguments
    data_raw Movie
    sigma double
end

fprintf('Applying Gaussian blur - sigma: %s\n',num2str(sigma))
    
% initialize
movie_result = data_raw;
movie = data_raw.stack;

for i = 1:data_raw.nfr
    movie(:,:,i) = imgaussfilt(movie(:,:,i),sigma);
end

% store to results
movie_result.stack = movie;

end
