function A = MakeRef34MAT(stim_sequence,pl)

A = zeros(numel(stim_sequence));
stims = unique(stim_sequence);
odortrial_idx = [];
for i_stim = 1:numel(stims)
    odortrial_idx(i_stim,:) = ismember(stim_sequence,stims{i_stim});
end
A = odortrial_idx'*odortrial_idx;

%%
if pl
figure; imagesc(A)

sortorder = 1:35;
labs = stim_sequence(sortorder);

    axis square; hold on
    xticks(1:max(sortorder)); xticklabels(labs(sortorder)); xtickangle(90)
    yticks(1:max(sortorder)); yticklabels(labs(sortorder))
    % xlabel('Stimulus')
    % ylabel('Stimulus')
    colorbar
    colormap(['hot'])
    colorbar('Color','w')
    set(gca, 'color', 'none', 'XColor','w', 'YColor','w', 'ZColor','w');
    set(gcf, 'color', 'none'); 
    hold off
end
end