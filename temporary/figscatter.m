function [h,FileOut] = figscatter(y,figlabs,do_animate)

global ntrials X labs pathout C trialn2use s stims2use stims
% scatterplot

%colormap
if numel(stims2use) == 1
    c = parula;
    c = c(floor(linspace(1,size(c,1),ntrials)),:);
else
    c = C(X,:);
end

scalefact = 500;
sizesc = 10+scalefact/(size(y,1)*size(y,3));

%plot
h = figure;
for i_trial = 1:size(y,3)
    switch size(y,2)
        case 2
            hs = scatter(y(:,1,i_trial),y(:,2,i_trial), ...
                sizesc,c(i_trial,:),'filled');
        case 3
            hs = scatter3(y(:,1,i_trial),y(:,2,i_trial),y(:,3,i_trial), ...
                sizesc,c(i_trial,:),'filled');
        otherwise
            error('dimension number not supported')
    end
    if numel(stims2use) == 1
        hs.DisplayName = ['trial #', num2str(trialn2use(i_trial))];
    else
        hs.DisplayName = stims{X(i_trial)};
    end
    hold on
end

% cosmetics / labels
title(figlabs{1},'color','w')
u = legendUnq();
legend(u,'Box','on','color','none','Location','best','EdgeColor','w','TextColor','w')
set(gcf, 'color', 'none');    
set(gca, 'color', 'none', 'XColor','w', 'YColor','w', 'ZColor','w');
xlabel(figlabs{2})
ylabel(figlabs{3})
if size(y,2)>2; zlabel(figlabs{4}); end
grid off
axis square

% save
% FileOut = fullfile(pathout,strcat('UMAPscatter','.fig'));
% if s; savefig(h,FileOut); end
% FileOut = fullfile(pathout,strcat(figlabs{1},'.gif'));
% if size(y,2)>2 && do_animate; animate3D(FileOut,3); end
rotate3d on