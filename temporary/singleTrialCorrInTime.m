function data = singleTrialCorrInTime(data, stim_on, stim_off, pl, s)
% singleTrialCorrInTime(data, stim_on, stim_off, pl, s)
% stim_on and stim_off in seconds


%% corr matrix between istantaneous activity vectors across timeframes

% input example: data, stim_on, stim_off
% indata= data.traces;
indata = data.tracesdns;
ncells = data.N;
indur = data.Ldns; % single trial duration in frames (after filtering)
intrialnums = data.trial_num;
intrialname = data.stim_type;
fs = data.meta.framerate/data.meta.downsample;
trials = data.trials;

ntrials = numel(trials)


figlab = {'cosinedistance'};

for i_fig = 1:numel(figlab)
    thisfiglab = figlab{i_fig};
    
    for i_trial = 1:ntrials
        %% Execute
        ttlstr = ['Trial #',num2str(intrialnums(i_trial)),' - ',intrialname{i_trial}];
        fprintf([ttlstr,' ... '])
        
        thistrialtraces = reshape(indata(:,i_trial), indur, ncells);
        %figure; imagesc(thistrialtraces)
        
        switch i_fig
            case 1
                thismetric = 'cosine';
                do_pca = true;
        end
        
        c= squareform(pdist(thistrialtraces,thismetric));
        
        
        if pl; hf = {}; end
        t = 0:indur/fs+1/fs;
        
        
        if do_pca
            Mpca = pca([thistrialtraces(2:end,:)]'); 
            if pl
                hf{end+1} = figure; imagesc(1:10,t(2:end),Mpca(:,1:min(size(Mpca,2),10))); hold on
                xlabel('PC#'); ylabel('time [s]')
                title([ttlstr,' - PCA (1st 10 PCs)'])
            end
        end
        
        if pl
            hf{end+1} = figure; imagesc(t,t,c), axis square, colorbar, hold on
            caxis([min(caxis), quantile(c(:), .999)])
            t = stim_on; line([0,t,t], [t,t,0], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1)
            t = stim_off; line([0,t,t], [t,t,0], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1)
            xlabel('time [s]'), ylabel('time [s]')
            title([ttlstr,' - ', figlab{i_fig}])
        end
        
        diag = indur * ([1:indur] - 1) + [1:indur]; % diagonal indices
        
        dxidt = c(diag(1:end-1)+1);
        dxidt_smooth = movmean(dxidt,3);
        t = linspace(0,indur/fs+1/fs, length(dxidt_smooth));
        if pl; hf{end+1} = figure; plot(t,dxidt_smooth); end
        
        
        c1 = c;
%         c1(diag(1:end-1)) = c1(diag(1:end-1)+1);
%         c1(end) = c1(end-1);
        intv_start = stim_on-5;
        intv_end = stim_off+10;
        interval = floor(intv_start*fs):floor((intv_end)*fs);
        c1 = c1(interval, interval); % stimulus window only
%         [clust,idx] = sort(clusterdata(c1,6));

        clustmethod = 'adaptAPclust';
        [labels, NCopt, NCopt2, clustparam] = adaptAPclust(thistrialtraces(interval,:));
        doclust = true;
        if NCopt==0; warning('Clustering didn''t work. Skipped.'); doclust = false; end
        if NCopt2==0; NCopt2 = NCopt; end
        
        if doclust
             %% correct cluster assignment

             % optimal solution
            [a, clustlab] = correctClustAssignment(labels, NCopt)
            nclustcorrect = numel(clustlab);
            labelscorrect = a';

            % second optimal solution
            [a, clustlab] = correctClustAssignment(labels, NCopt2)
            nclustcorrect2 = numel(clustlab);
            labelscorrect2 = a';


            %% plot clustering

            if pl

                % in chronological order
                t = intv_start - stim_on + [0:size(c1,1)/fs+1/fs];
                y = 1-c1;
                hf{end+1} = figure; imagesc(t,t,y), axis square, colorbar, hold on
                caxis([min(caxis), quantile(y(:), .999)])
                title([ttlstr ' - interval of interest'])
                xlabel('time from stimulus onset [s]')
                ylabel('time from stimulus onset [s]')

                % by cluster opt
                T = labels(:,min(size(labels,2),NCopt-1)); [clust,idx] = sort(T);
                hf = plotintervalbyclust(hf, clust, idx, c1, t);
                title([ttlstr ' - ap clust'])

                % by cluster opt 2 
                T = labels(:,min(size(labels,2),NCopt2-1)); [clust,idx] = sort(T);
                hf = plotintervalbyclust(hf, clust, idx, c1, t);
                title([ttlstr ' - ap clust2'])

                % by corrected cluster opt 1 
                clust = labelscorrect; idx = 1:numel(labelscorrect);
                hf = plotintervalbyclust(hf, clust, idx, c1, t);
                title([ttlstr ' - ap clust corrected'])

                % by corrected cluster opt 1 
                clust = labelscorrect2; idx = 1:numel(labelscorrect2);
                hf = plotintervalbyclust(hf, clust, idx, c1, t);
                title([ttlstr ' - ap clust2 corrected'])

            end

            %%
        end
        
        
        %% saving figures
        if pl && s
            fprintf('Saving ... ')
            
            if ~exist(mfilename,'dir'); mkdir(mfilename); end
            
            nameflag = [replace(data.stim_type{i_trial},'.',''),'_',num2str(data.trial_num(i_trial)),'_'];
            
            fnames = {'pca1to10', thisfiglab, 'dxidtsmooth', 'ioi', ...
                'apclust1', 'apclust2', 'apclustcorrect', 'apclustcorrect2'};
            
            for i = 1:numel(hf)
                pathout = [mfilename,'\',fnames{i}];
                if ~exist(pathout,'dir'); mkdir(pathout); end
                
                FileOut = [pathout,'\', nameflag, fnames{i}];
                %print([FileOut,'.pdf'],'-dpdf','-fillpage')
                savefig(hf{i},FileOut)
                saveas(hf{i},FileOut,'jpeg')
                
                close(hf{i})
            end
            
        end
        
        
        fprintf(' \n')
                     
        
        %% updating data (not saved to file but returned as var)
        
        out.pca = Mpca;
        out.distancemetric{i_fig} = thismetric;
        out.distancemat{i_fig} = c;
        out.ddistancedt = dxidt;
        
        out.resp_stage_clust{1}.window_sec = [intv_start; intv_end];
        out.resp_stage_clust{1}.methodname = clustmethod;
        out.resp_stage_clust{1}.methodparam = clustparam;
        out.resp_stage_clust{1}.clustlabels = labels;
        out.resp_stage_clust{1}.optNClust = NCopt;
        out.resp_stage_clust{1}.optNClust2 = NCopt2;
        out.resp_stage_clust{1}.nclustcorrect = nclustcorrect;
        out.resp_stage_clust{1}.labelscorrect = labelscorrect;
        out.resp_stage_clust{1}.nclustcorrect2 = nclustcorrect2;
        out.resp_stage_clust{1}.labelscorrect2 = labelscorrect2;
        
        
        
        data.singletrial{i_trial}.intime = out;
        

    end
    
    

end

fprintf('... DONE! \n')



end



function hf = plotintervalbyclust(hf, clust, idx, c1, t)

hf{end+1} = figure;
h = subplot(131); set(h,'Units','pixels','Position',[50 30 10 360])
imagesc(idx(:)); title('orig t')
h = subplot(132); set(h,'Units','pixels','Position',[100 30 10 360])
imagesc(clust(:)); title('clust')
h = subplot(133); set(h,'Units','pixels','Position',[150 30 360 360])
imagesc(t,t,c1(idx, :)), axis square, %colorbar, hold on
caxis([min(caxis), quantile(c1(:), .999)])

end


function [a, clustlab] = correctClustAssignment(labels, optNClust)


idx = find(max(labels)==optNClust);
a = [labels(:,idx)]';


% prune loner assignments
for i = 2:length(a)-1
    if a(i-1)~=a(i) && a(i+1)~=a(i)
        a(i) = a(i-1);
    end
end
a(1) = a(2);
a(end) = a(end-1);
a0 = a;

nperiods = numel(find(diff(a)~=0))+1;


% assign new names to split periods
if nperiods ~= optNClust
    disp(['Adjusting from ',num2str(optNClust),' to ', num2str(nperiods), ' stages.'])

    % initialize vars
    clustlab = [];
    convert = false;
    toconvert = 0;

    for i = 1:length(a)
        if convert
            if a(i)==toconvert; a(i) = a(i-1);
            else; convert = false;
            end
        end
        if i==1 || ( a(i-1)~=a(i) && ~ismember(a(i),clustlab) )
            clustlab = [clustlab; a(i)]; 
        elseif a(i-1)~=a(i) && ismember(a(i),clustlab)
            convert = true;
            newlab = 1;
            while ismember(newlab,clustlab)
                newlab = newlab+1; % identify first unused integer label
            end

            toconvert = a(i);
            a(i) = newlab;
            clustlab = [clustlab; a(i)]; 
        end
    end
else
    clustlab = unique(a);
end

clustlab

%  temp plot
%figure; imagesc([a', a0', labels'])

end
