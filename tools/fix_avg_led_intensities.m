function ai = fix_avg_led_intensities(p, f, overwrite, do_plot)
% Load, optionally skip if already corrected, then visualize "before",
% prompt for edits, interpolate, trim, save, and (optionally) show "after".
% p         folder path, e.g. 'tail_movies/rot/LED_vals'
% f         filename, e.g. 'r_avg_intensities.mat' (contains struct avg_intensities)
% overwrite true/false: if false and .allchecked exists & non-empty, return immediately
% do_plot   true/false to also show "after" plot (default true)

if nargin < 3 || isempty(overwrite), overwrite = false; end
if nargin < 4 || isempty(do_plot),   do_plot   = true;  end

s = load(fullfile(p, f), 'avg_intensities');
ai = s.avg_intensities;

% early exit if already corrected and not overwriting
if ~overwrite && isfield(ai, 'allchecked') && ~isempty(ai.allchecked)
    return
end

% show "before"
figure; plot(ai.all); title('Before: avg\_intensities.all'); xlabel('Frame'); ylabel('Avg intensity'); drawnow;

% prompt while viewing "before"
dlg = inputdlg( ...
    { 'Interpolation START index (empty to skip):', ...
      'Interpolation END index (empty to skip):', ...
      'Final frame to keep (empty to skip):' }, ...
    'Fix avg LED intensities', 1, {'', '', ''});

% parse inputs
r1 = []; r2 = []; endfr = [];
if ~isempty(dlg) % handle cancel
    if ~isempty(dlg{1}), r1 = str2double(dlg{1}); end
    if ~isempty(dlg{2}), r2 = str2double(dlg{2}); end
    if ~isempty(dlg{3}), endfr = str2double(dlg{3}); end
end

ai.allchecked = ai.all;

% apply interpolation window if valid
if ~isempty(r1) && ~isnan(r1) && ~isempty(r2) && ~isnan(r2)
    r1 = max(1, round(r1));
    r2 = min(numel(ai.allchecked), round(r2));
    if r2 >= r1
        ai.allchecked(r1:r2) = NaN;
    end
end

% fill NaNs linearly
ai.allchecked = fillmissing(ai.allchecked, 'linear');

% trim if requested
if ~isempty(endfr) && ~isnan(endfr)
    endfr = max(1, min(numel(ai.allchecked), round(endfr)));
    ai.allchecked = ai.allchecked(1:endfr);
end

% optional "after" view
if do_plot
    figure; plot(ai.allchecked); title('After: avg\_intensities.allchecked'); xlabel('Frame'); ylabel('Avg intensity');
end

% save back to file
s.avg_intensities = ai;
save(fullfile(p, f), '-struct', 's');
end
