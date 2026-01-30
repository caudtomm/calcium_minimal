function out = characterizePCspace(a, coef,explained, cfg, nshuffle)
arguments
    a double
    coef double
    explained double
    cfg PlotConfig = PlotConfig()
    nshuffle double = 50;
end

% how much of the variance does each PC explain?
shuffleexplained = [];
for i_sh = 1:nshuffle
    shufflecoef = coef(randperm(size(coef,1)),:);
    shufflescore = a*shufflecoef;
    shufflepcvar = std(shufflescore,[],1,'omitnan').^2;
    shuffleexplained(i_sh,:) = (shufflepcvar*100)./sum(shufflepcvar);
end
shufflemean = nanmean(shuffleexplained,1);
shufflestd = std(shuffleexplained);
pcthr = movmean(shufflemean + shufflestd*2,5)';
out.maxPC = find(explained<=pcthr,1)-1;

out.hf = figure; hold on
y = shufflemean; t = 1:length(y);
curve1 = y + shufflestd; curve2 = y - shufflestd;
h = patch([t,fliplr(t)],[curve1, fliplr(curve2)],'g','FaceAlpha',.3,'EdgeColor','none');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot(t,y,'Color','g','LineWidth',2)
plot(t,explained,'Color','b','LineWidth',2)
plot(t,pcthr,'Color','r','LineStyle','--','LineWidth',1)
xlim([t(1),t(end)])
axis square
set(gcf, 'color', cfg.bgcol);    
set(gca, 'color', cfg.bgcol', 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
legend({'shuffle','data','thres (2 STD)'},'Box','on','color',cfg.bgcol,'Location','best','EdgeColor',cfg.axcol,'TextColor',cfg.textcol)
ylabel('Variance explained')
xlabel('PC #')

end