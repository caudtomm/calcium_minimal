
figure; imagesc(mean(singlenorm,3))

title('single trial norm(qt .03)-avg 3x OB')
xticks([1:11])
xticklabels({'Ser','Ser','Ser', ...
    'Ser+Var(15'')','Ser+Var(30'')','Ser+Var(60'')','Ser+Var(90'')',...
    'Ser(15'')','Ser(30'')','Ser(60'')','Ser(90'')',})
xtickangle(45)
t0 = linspace(0,311,8)
yticks(t0(2:7))
t = linspace(25,65,8)
yticklabels(t)
yticklabels(t(:))
yticklabels(floor(t(2:7)))
line([0,12],[30*7.67-189,30*7.67-189],'LineWidth',2,'Color','red','LineStyle','--')
line([0,12],[50*7.67-189,50*7.67-189],'LineWidth',2,'Color','red','LineStyle','--')
ylabel('time [s]')
set(gcf, 'Position', [50 50 250 700]);


figure; imagesc(mean(globalnorm,3))

title('global trial norm(qt .03)-avg 3x OB')
xticks([1:11])
xticklabels({'Ser','Ser','Ser', ...
    'Ser+Var(15'')','Ser+Var(30'')','Ser+Var(60'')','Ser+Var(90'')',...
    'Ser(15'')','Ser(30'')','Ser(60'')','Ser(90'')',})
xtickangle(45)
t0 = linspace(0,311,8)
yticks(t0(2:7))
t = linspace(25,65,8)
yticklabels(t)
yticklabels(t(:))
yticklabels(floor(t(2:7)))
line([0,12],[30*7.67-189,30*7.67-189],'LineWidth',2,'Color','red','LineStyle','--')
line([0,12],[50*7.67-189,50*7.67-189],'LineWidth',2,'Color','red','LineStyle','--')
ylabel('time [s]')
set(gcf, 'Position', [50 50 250 700]);