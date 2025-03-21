%%

pl = 0;

% get response intensity for each cell and trial
interval_base = [floor((data.f0_window(1))*fs) : (data.f0_window(2))*fs];
interval_response = [floor((data.stim_on_sec+1)*fs) : (data.stim_off_sec)*fs];
thistraces = traceFormat(data.traces,data.L);
cellstd_base = mean(squeeze(std(thistraces(interval_base,:,:),[],1,'omitnan')),2,'omitmissing');
cellresponse = squeeze(mean(thistraces(interval_response,:,:),1,'omitmissing'));% ./ cellstd_base;


% get variance for each cell and trial
interval = [floor((data.stim_on_sec-5)*fs) : (data.stim_off_sec)*fs];
thistraces = traceFormat(data.traces,data.L);
cellstd = squeeze(std(thistraces(interval,:,:),[],1,'omitnan'));
cellvar = cellstd.^2;


if pl
    %------------
    figure;
    subplot(121), imagesc(cellvar)
    title('cell variance')
    xticks(1:size(cellvar,2)); xticklabels(data.stim_type)
    ylabel('cell #')
    subplot(122), imagesc(cellvar(:,data.idx_by_stim_type))
    title('cell variance')
    xticks(1:size(cellvar,2)); xticklabels(data.stim_type(data.idx_by_stim_type))
    %------------
    figure
    violin(cellvar)
    xticks(1:34); xticklabels(data.stim_type); xtickangle(90)
    ylabel('Unit variance during stimulus')
    %------------
    figure
    violin(cellvar(:,data.idx_by_stim_type))
    xticks(1:34); xticklabels(data.stim_type(data.idx_by_stim_type)); xtickangle(90)
    ylabel('Unit variance during stimulus')
    %------------
end




[h,idx] = sort(cellvar); %idx = idx(end-19:end,:);

if pl
%------------
figure; loglog(h./max(h))
axis tight, box off
xlabel('unit ID (sorted)')
ylabel('normalized unit variance')
%------------
figure; loglog(mean(h))
axis tight, box off
ylabel('mean unit variance')
xlabel('trial number')
%------------
figure; plot(cumsum(h)./max(cumsum(h)))
axis tight, box off
xlabel('unit ID (sorted)')
ylabel('normalized cumulative variance')
%------------
figure; loglog(cumsum(h)./max(cumsum(h)))
axis tight, box off
xlabel('unit ID (sorted)')
ylabel('normalized cumulative variance')
%------------
end



P = .8;

% get units contributing to P% of variance
ncsum = cumsum(h)./max(cumsum(h));
topunits = {}; ntopunits = [];
for i_trial = 1:size(ncsum,2)
    thr = find(ncsum(:,i_trial)>P, 1); % 80% of variance
    topunits{end+1} = idx(thr:end,i_trial);
    ntopunits = [ntopunits numel(topunits{end})];
end

if pl
%------------
figure; plot(ntopunits/data.N,'LineWidth',2)
xticks(1:34); xticklabels(data.stim_type); xtickangle(90)
ylabel('portion of top 80% variant units')
box off, ylim([0,1])
%------------
end

topunits2 = topunits(data.idx_by_stim_type);
M = []; M2 = [];
for i_trial = 1:size(ncsum,2)
    for i_trial2 = 1:size(ncsum,2)
        thisintersect = intersect(topunits{i_trial},topunits{i_trial2});
        thisintersect2 = intersect(topunits2{i_trial},topunits2{i_trial2});
        M(i_trial,i_trial2) = numel(thisintersect)/numel(topunits{i_trial});
        M2(i_trial,i_trial2) = numel(thisintersect2)/numel(topunits2{i_trial});
    end
end

if pl
%------------
figure; imagesc(M), axis square, colorbar
xticks(1:34); xticklabels(data.stim_type); xtickangle(90)
yticks(1:34); yticklabels(data.stim_type); xtickangle(90)
title('portion of common top variant units')
%------------
figure; imagesc(M2), axis square, colorbar
xticks(1:35); xticklabels(data.stim_type(data.idx_by_stim_type)); xtickangle(90)
yticks(1:35); yticklabels(data.stim_type(data.idx_by_stim_type)); xtickangle(90)
title('portion of common top variant units')
%------------
end


%%

alltopunits = [];
for i=1:numel(topunits)
    alltopunits = [alltopunits;topunits{i}];
end
alltopunits = unique(alltopunits);

