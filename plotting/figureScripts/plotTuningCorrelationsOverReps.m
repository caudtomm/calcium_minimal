function [h, C, Cvals] = plotTuningCorrelationsOverReps(v, metric)

cfg = v.plotConfig;

tc = v.plotUnitActivityMetricHead('plotType', 'none', ...
                        'method', 'tuning curves', ...
                        'n_equals', 'cells', ...
                        'do_normalize', false); % cell, each [units x stims x reps]

tc = cell2mat(tc);
[nUnits, nstims, nReps] = size(tc);

C = zeros(nReps, nReps,nUnits);
Cvals = zeros(nUnits,1);
for i = 1:nUnits
    thisunit = squeeze(tc(i,:,:)); % [nStims x nReps]
    C(:,:,i) = 1 - squareform(pdist(thisunit',metric));
    Cvals(i) = 1 - mean(pdist(thisunit',metric),'all','omitmissing');
end

h = figure;

subplot(131)
histogram(Cvals,'FaceColor','k')
hold on
line([.5,.5],[0 300],'Color','r','LineStyle','--','LineWidth',2)
set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
axis square
box off
xlabel(['avg inter-rep ',metric])

subplot(132)
imagesc(mean(C,3,'omitmissing')); axis square
a = colorbar('Color',cfg.axcol);
a.Label.String = metric;
a.Label.FontSize= gca().FontSize;
title('all units or modes')

subplot(133)
idx = Cvals>.5;
imagesc(mean(C(:,:,idx),3,'omitmissing')); axis square
a = colorbar('Color',cfg.axcol);
a.Label.String = metric;
a.Label.FontSize= gca().FontSize;
title(['avg inter-rep',metric,' >50%'])

set(gcf, 'color', cfg.bgcol);
set(gcf,'position',[1 1 1000 500])