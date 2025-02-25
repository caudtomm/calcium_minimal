function [h,FileOut] = fig3curves(y,figlabs,do_patch)
% y : [t x dim x trial]
% figlabs : {'title', 'label1', 'label2', 'label3'}
global ntrials X labs C pathout s s2u_idx tframe fs trialn2use

if ~exist("do_patch","var") || isempty(do_patch); do_patch = false; end

h = figure;
subplot(311); hold on; ylabel(figlabs{2},'Color','w')
subplot(312); hold on; ylabel(figlabs{3},'Color','w')
subplot(313); hold on; ylabel(figlabs{4},'Color','w'); xlabel('time from onset [s]','Color','w')

t = tframe(1) + [0:1/fs:size(y,1)/fs-1/fs];

c = parula;
c = c(floor(linspace(1,size(c,1),ntrials)),:);

if ~do_patch

    for i_trial = 1:size(y,3)
        if numel(unique(X)) == 1
            thisc = c(i_trial,:);
            thisname = ['trial #', num2str(trialn2use(i_trial))];
        else
            thisc = C(X(i_trial),:);
            thisname = labs{X(i_trial)}';
        end
        subplot(311); plot(t,y(:,1,i_trial),'Color',thisc,'DisplayName',thisname);
        subplot(312); plot(t,y(:,2,i_trial),'Color',thisc);
        subplot(313); plot(t,y(:,3,i_trial),'Color',thisc);
    end

else

    for i_stim = 1:numel(s2u_idx)
        idx = ismember(X,s2u_idx(i_stim));
        tmp = y(:,:,idx);
        stimmean = squeeze(nanmean(tmp,3));
        stimstd = squeeze(std(tmp,[],3,'omitnan'));
    
        ymean = stimmean;
        curve1 = ymean + stimstd; curve2 = ymean - stimstd;
        
        for i_dim = 1:3
            subplot(3,1,i_dim); 
            hp = patch([t,fliplr(t)],[curve1(:,i_dim); fliplr(curve2(:,i_dim)')'], ...
                C(s2u_idx(i_stim),:),'FaceAlpha',.3,'EdgeColor','none');
            hp.Annotation.LegendInformation.IconDisplayStyle = 'off';
            plot(t,ymean(:,i_dim),'Color',C(s2u_idx(i_stim),:),'LineWidth',2,'DisplayName',labs{s2u_idx(i_stim)});
        end
    end

end

for i_sp = 1:3
    subplot(3,1,i_sp); 
    if i_sp == 1
        u = legendUnq();
        legend(u,'Box','off','color','none','Location','best','EdgeColor','w','TextColor','w')
    end
    set(gca, 'color', 'none', 'XColor','w', 'YColor','w', 'ZColor','w');
    grid off
    axis tight
end
set(gcf, 'color', 'none'); 
set(gcf, 'Position', [50 50 700 800]);

FileOut = fullfile(pathout,strcat(figlabs{1},'.fig'));
if s; savefig(h,FileOut); end
FileOut = fullfile(pathout,strcat(figlabs{1},'.svg'));
if s; plot2svg(FileOut,h); end