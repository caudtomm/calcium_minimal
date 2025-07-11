function windows = defineTimeWindows(window_duration, t_lim_sec, overlap_sec)
%DEFINE_TIME_WINDOWS Generate overlapping time windows in seconds.
%
%   WINDOWS = DEFINE_TIME_WINDOWS(WINDOW_DURATION, T_LIM_SEC, OVERLAP_SEC)
%   returns an N-by-2 array WINDOWS where each row specifies the [START END]
%   time (in seconds) of a window. Windows span the time range defined by
%   T_LIM_SEC and are WINDOW_DURATION seconds long, overlapping by
%   OVERLAP_SEC seconds.
%
%   Inputs:
%     - WINDOW_DURATION : duration of each window in seconds (scalar > 0)
%     - T_LIM_SEC       : 1x2 vector [start_time end_time] in seconds
%     - OVERLAP_SEC     : overlap between consecutive windows in seconds
%                         (scalar >= 0 and < WINDOW_DURATION)
%
%   Output:
%     - WINDOWS         : Nx2 matrix, each row is [start_time end_time] in sec

    % Validate inputs
    if nargin < 3
        overlap_sec = 0;
    end
    if nargin < 2
        t_lim_sec = [-1, 25];
    end
    if nargin < 1
        window_duration = 1;
    end

    assert(window_duration > 0, 'window_duration must be positive.');
    assert(overlap_sec >= 0 && overlap_sec < window_duration, ...
        'overlap_sec must be non-negative and smaller than window_duration.');
    assert(numel(t_lim_sec) == 2 && t_lim_sec(1) < t_lim_sec(2), ...
        't_lim_sec must be a 2-element vector [start end] with start < end.');

    % Step size is reduced by overlap
    step = window_duration - overlap_sec;

    % Generate window start times
    t_start = t_lim_sec(1);
    t_end = t_lim_sec(2);
    starts = t_start:step:(t_end - window_duration);
    
    % Construct windows
    windows = [starts(:), starts(:) + window_duration];
end
