% example: first Trp trial

global tframe

stims2use = {'Arg','Ala','His','Trp','Ser','Leu','ACSF','spont.'};
nstims = numel(stims2use);
trialn2usein = 1:12;
s = 0;
plot_intermediate = false;

todo_fish = 1;
todo_trials = 1:30;

method = 'cosine';
mode = 'forward'; % {'forward', 'reverse', 'insensitive'}


tframe = [0 0];

%%
[a, X,stims,~,~, trialidx] = selectData(experiment, todo_fish, stims2use, trialn2usein,'trialavg');
ntrials = numel(X);

%% make optimal corr mat
A = zeros(35);
stims = unique(X);
odortrial_idx = [];
for i_stim = 1:numel(stims)
    odortrial_idx(i_stim,:) = ismember(X,stims(i_stim));
end
A = odortrial_idx'*odortrial_idx;

% A = X'*X;
sortorder = trialidx;
labs = experiment.series{todo_fish(end)}.data.stim_type(sortorder);

c = A;
figure; imagesc(c)
axis square; hold on
xticks(1:size(c,1)); xticklabels(labs); xtickangle(90)
yticks(1:size(c,1)); yticklabels(labs)
xlabel('Trial number')
ylabel('Stimulus type')
colorbar
hold off

%% action

% stims = {'spont.','ACSF','Trp','Ser','Ala','Food'};
stims = {'Arg','Ala','His','Trp','Ser','Leu'};
niter = 10;
[ntrials, ncells] = size(a);
kradius = 5;
nksamples = 100;
K = kradius./([1:floor(nksamples/2)].^2); K = [-K,fliplr(K)]';
learn_rate = .2;

K_opt = zeros(ncells,niter);

atmp = a;
atmpMAT = nan(ntrials,ncells,niter+1);
atmpMAT(:,:,1) = a;

y = nan(numel(todo_trials),ncells);

for i_trial = 1:numel(todo_trials)
    switch mode
        case {'forward', 'insensitive'}
            thistrial = todo_trials(i_trial);
        case 'reverse'
            thistrial = todo_trials(end-i_trial+1);
        otherwise
            error('specified mode is unknown')
    end
    disp(''); disp(['new trial: #',num2str(thistrial)]); disp('')
    for i_iter = 1:niter
        disp(['iteration #',num2str(i_iter),' of ',num2str(niter)])
        for i_cell = 1:ncells
            Distancebyk = [];
            for j_k = 1:nksamples
                if ~strcmp(mode,'insensitive')
                    avg_actvect = a;
                else
                    avg_actvect = atmp;
                end
                avg_actvect(thistrial,i_cell) = avg_actvect(thistrial,i_cell) * (1 + K(j_k));
                thiscorrmat = 1-squareform(pdist(avg_actvect, method));
                thisdistance = nanmean(abs(A(thistrial,:)-thiscorrmat(thistrial,:)));
                Distancebyk = [Distancebyk; thisdistance];
            end
            [~, idx] = min(Distancebyk); 
            K_opt(i_cell,i_iter) = K(idx);
            atmp(thistrial,i_cell) = atmp(thistrial,i_cell) .* (1 + K_opt(i_cell,i_iter) * learn_rate)';
        end
        
        atmpMAT(:,:,i_iter+1) = atmp;
    end
    
    if plot_intermediate
        figure; 
        for i=1:niter+1
            thiscorrmat = 1-squareform(pdist(atmpMAT(:,:,i), method));
            imagesc(thiscorrmat)
            caxis([0 1])
            yticks(1:ntrials); yticklabels(stims(X))
            xticks(1:ntrials); xticklabels(stims(X))
            colormap('copper')
            colorbar('Color','w')
            set(gca, 'color', 'none', 'XColor','w', 'YColor','w', 'ZColor','w');
            set(gcf, 'color', 'none'); 
            title(num2str(i-1),'color','w')
            axis square
            pause(.5)
        end
    end

    y(i_trial,:) = sum(K_opt,2).*a(thistrial,:)'*learn_rate;
    if plot_intermediate
        figure; histogram(y(i_trial,:),100,'FaceColor','r')
        set(gca, 'color', 'none', 'XColor','w', 'YColor','w', 'ZColor','w');
        set(gcf, 'color', 'none'); box off
        xlabel('relative intensity change')
        ylabel('histogram')
        set(gcf, 'Position', [50 50 350 200]);
    end
end

figure; histogram(y,100,'FaceColor','r')
set(gca, 'color', 'none', 'XColor','w', 'YColor','w', 'ZColor','w');
set(gcf, 'color', 'none'); box off
xlabel('relative intensity change')
ylabel('histogram')
set(gcf, 'Position', [50 50 350 200]);

%% plotting

figure; imagesc(y)
caxis([-1 1])
z = std(y,[],1).^2;
z_idx = z>quantile(z,.8);
figure; imagesc(z_idx)
figure; imagesc(exp(y(:,z_idx)))
caxis([-1 1])


%%

[Xsort,idx] = sort(X);

Xsort = X; idx = 1:length(Xsort);

figure; imagesc(y(idx,:)); caxis([-1 1])
yticks(1:ntrials); yticklabels(stims(Xsort))

changestd = std(y(idx,:),[],1);
[~,changesort] = sort(changestd);
figure; imagesc(y(idx,changesort)); caxis([-1 1])
yticks(1:ntrials); yticklabels(stims(Xsort))

figure; imagesc(zscore(y(idx,changesort)')'); caxis([-4 4])
yticks(1:ntrials); yticklabels(stims(Xsort))

figure
thiscorrmat = 1-squareform(pdist(y(Xsort,:), method));
imagesc(thiscorrmat)
yticks(1:ntrials); yticklabels(stims(Xsort))
xticks(1:ntrials); xticklabels(stims(Xsort))
caxis([0 1])
colormap('hot')
colorbar('Color','w')
set(gca, 'color', 'none', 'XColor','w', 'YColor','w', 'ZColor','w');
set(gcf, 'color', 'none'); 
axis square





