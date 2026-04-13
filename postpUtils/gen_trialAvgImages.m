function gen_trialAvgImages(movie_dir, frame_start, frame_end, out_file)
% gen_trialAvgImages  Compute per-trial mean dF/F images over a frame interval.
%
%   gen_trialAvgImages(movie_dir, frame_start, frame_end, out_file)
%
%   Loads all TC_*.mat dF/F movie files from movie_dir in sorted order,
%   averages each over frame_start:frame_end (omitting NaN frames), and
%   stacks the results into a [H x W x nTrials] double saved to out_file.
%
%   movie_dir    directory containing per-trial dF/F .mat files (TC_*.mat)
%   frame_start  first frame of the averaging interval (1-based, inclusive)
%   frame_end    last frame of the averaging interval (inclusive)
%   out_file     output .mat file path; variable name in file is 'stack'

files = dir(fullfile(movie_dir, 'TC_*.mat'));
files = files(~[files.isdir]);
[~, order] = sort({files.name});
files = files(order);
nfiles = numel(files);

if nfiles == 0
    error('gen_trialAvgImages:noFiles', 'No TC_*.mat files found in:\n  %s', movie_dir);
end

fprintf('Source  : %s  |  %d files\n', movie_dir, nfiles);
fprintf('Frames  : %d - %d\n', frame_start, frame_end);

stack = [];

for k = 1 : nfiles
    in_file = fullfile(files(k).folder, files(k).name);
    fprintf('[%d/%d] %s ... ', k, nfiles, files(k).name);

    movie    = robust_io('load', in_file, 'movie').movie;
    fr_range = frame_start : frame_end;
    fr_range = fr_range(fr_range >= 1 & fr_range <= movie.nfr);

    if isempty(fr_range)
        warning('gen_trialAvgImages:intervalOutOfRange', ...
            'Interval [%d %d] out of range for %s (nfr=%d) — inserting NaN.', ...
            frame_start, frame_end, files(k).name, movie.nfr);
        img = nan(movie.h, movie.w);
    else
        img = mean(movie.stack(:, :, fr_range), 3, 'omitmissing');
    end

    stack = cat(3, stack, img);
    fprintf('done\n');
end

out_dir = fileparts(out_file);
if ~isempty(out_dir) && ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

fprintf('Saving [%d x %d x %d] to: %s\n', size(stack,1), size(stack,2), size(stack,3), out_file);
save(out_file, 'stack', '-v7.3');
fprintf('Done.\n');
