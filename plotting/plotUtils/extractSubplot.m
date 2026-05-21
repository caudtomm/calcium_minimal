function hf_new = extractSubplot(hf, subplot_id)
% Return a new figure containing only subplot SUBPLOT_ID from figure HF.
% SUBPLOT_ID is 1-based in visual panel order (top-left first).

ax_all = flipud(findobj(hf, 'Type', 'axes'));
ax     = ax_all(subplot_id);

hf_new  = figure('Color', hf.Color);
new_ax  = copyobj(ax, hf_new);
new_ax.Position = [0.13, 0.11, 0.775, 0.815];

if ~isempty(ax.Legend)
    copies = copyobj([ax.Legend, ax], hf_new);
    delete(copies(2));  % remove duplicate axes; keep only the legend copy
end
end
