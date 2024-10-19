function idx_out = sortbyStimType(data)

a = unique(data.stim_type);

idx_4 = [];
for i_stim = 1:numel(a)
    idx = find(strcmp(data.stim_type,a{i_stim}));
    [~, idx_2] = sort(data.trial_num(idx));
    idx_3 = idx(idx_2);
    idx_4 = [idx_4; idx_3];
end

% output
idx_out = idx_4;