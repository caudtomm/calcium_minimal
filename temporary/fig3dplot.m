function [h,FileOut] = fig3dplot(y,figlabs,linew,do_animate)
% y : [t x dim x trial]
% figlabs : {'title', 'xlabel', 'ylabel', 'zlabel'}
% (linew) : double, numel=ntrials
global ntrials X labs pathout C trialn2use s
if ~exist("linew",'var') || isempty(linew); linew = repelem(2,ntrials); end
h = figure;
c = parula;
c = c(floor(linspace(1,size(c,1),ntrials)),:);
thisstims = unique(X);
for i_trial = 1:ntrials
    if size(y,2) == 2
        hp = plot(y(:,1,i_trial),y(:,2,i_trial), ...
            'Color',C(X(i_trial),:), ...
            'LineWidth',linew(i_trial), ...
            'DisplayName',labs{X(i_trial)});
    elseif size(y,2) == 3
        hp = plot3(y(:,1,i_trial),y(:,2,i_trial),y(:,3,i_trial), ...
            'Color',C(X(i_trial),:), ...
            'LineWidth',linew(i_trial), ...
            'DisplayName',labs{X(i_trial)});
    else
        error('number of dimensions not supported for plotting')
    end

    if numel(thisstims) == 1
        hp.Color = c(i_trial,:);
        hp.DisplayName = ['trial #', num2str(trialn2use(i_trial))];
    end
    hold on
end
% for i_stim = 1:numel(thisstims)
%     idx = ismember(X,thisstims(i_stim));
%     ymean = nanmean(y(:,:,idx),3);
%     hp = plot3(ymean(:,1),ymean(:,2),ymean(:,3), ...
%             'Color',C(thisstims(i_stim),:), ...
%             'LineWidth',2, ...
%             'DisplayName',[labs{thisstims(i_stim)}, ' - AVG']);
%     hold on
% end
title(figlabs{1})
u = legendUnq();
legend(u,'Box','on','color','none','Location','best','EdgeColor','w','TextColor','w')
xlim([min(xlim)-1,max(xlim)+1])
ylim([min(ylim)-1,max(ylim)+1])
zlim([min(zlim)-1,max(zlim)+1])
set(gcf, 'color', 'none');    
set(gca, 'color', 'none', 'XColor','w', 'YColor','w', 'ZColor','w');
grid on
axis square
xlabel(figlabs{2})
ylabel(figlabs{3})
zlabel(figlabs{4})
FileOut = fullfile(pathout,strcat(figlabs{1},'.fig'));
if s; savefig(h,FileOut); end
FileOut = fullfile(pathout,strcat(figlabs{1},'.gif'));
if do_animate; animate3D(FileOut,3); end
rotate3d on

