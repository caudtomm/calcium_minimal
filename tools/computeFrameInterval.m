function frame_interval = computeFrameInterval(traces, trial_num, interval_sec)
% computeFrameInterval  Convert a [t_start, t_end] interval (seconds relative
%   to effective stimulus onset) into movie frame indices.
%
%   frame_interval = computeFrameInterval(traces, trial_num, interval_sec)
%
%   Effective stimulus onset = stim_series.frame_onset + odor_delay * fs.
%   trial_num is clamped to the valid row range of stim_series.
%
%   traces        ActivityTraces
%   trial_num     scalar — 1-based trial index (row into stim_series)
%   interval_sec  [t_start, t_end] in seconds relative to stimulus onset

fs            = traces.framerate;
row           = min(trial_num, height(traces.stim_series));
delay         = traces.odor_delay; if isempty(delay); delay = 0; end
stim_onset_fr = traces.stim_series.frame_onset(row) + round(delay * fs);
fr_start      = stim_onset_fr + round(interval_sec(1) * fs);
fr_end        = stim_onset_fr + round(interval_sec(2) * fs);
frame_interval = fr_start : fr_end;
