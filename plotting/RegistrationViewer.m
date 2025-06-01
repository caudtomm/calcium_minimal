classdef RegistrationViewer
    properties
        sj
        folder
    end

    methods
        function obj = RegistrationViewer(subject, foldername)
            arguments
                subject Subject
                foldername char = pwd
            end

            obj.sj = subject;            
            obj.folder = foldername;
        end

        function collage = subsampledMovieCollage(obj)
            orig_folder = fullfiletol(pwd);

            % get file list and move to subject-datapath
            [fileList, numFiles] = obj.sj.getFileListInSubfolder(obj.folder);
        
            % Calculate the dimensions of the collage
            cols = ceil(sqrt(numFiles));
            rows = ceil(numFiles / cols);
        
            % Load the first file to get the size of movie.stack
            disp('Initializing...')
            movie = robust_io('load',fileList{1}, 'movie').movie;
            [x, y, z] = size(movie.stack);
            
            % Subsample dimensions
            x = ceil(x / 5);
            y = ceil(y / 5);
        
            % Initialize the collage matrix
            collage = zeros(rows * x, cols * y, z);
            
            % Add each movie to the collage
            for i = 1:numFiles
                % Load the movie data
                disp(fileList{i})
                movie = robust_io('load',fileList{i}, 'movie').movie;
                
                % Subsample the movie stack
                subsampledStack = movie.stack(1:5:end, 1:5:end, :);
                
                % Determine the position in the collage
                rowIdx = ceil(i / cols);
                colIdx = mod(i-1, cols) + 1;
                
                % Insert the subsampled movie stack into the collage
                collage((rowIdx-1)*x + (1:x), (colIdx-1)*y + (1:y), :) = subsampledStack;
            end

            collage = Movie(collage);
            collage = collage.setFramerateHz(obj.sj.framerate);
            collage.path.fname = 'sidebyside';
            
            % return to original location
            cd(orig_folder)
        end

        function replaceFailsInterp(obj,fail_idx, sourceFolder, pl)
            % fail_idx : putative_fails output of QC()

            orig_folder = fullfiletol(pwd);

            % get file list and move to subject-datapath
            [fileList, numFiles] = obj.sj.getFileListInSubfolder(obj.folder);
            [fileList2, numFiles2] = obj.sj.getFileListInSubfolder(sourceFolder);
            if numFiles~=obj.sj.getNTrials; warning( ...
                    'Number of movies found  in the destination folder is different from the number of trials.'); end
            if numFiles2~=obj.sj.getNTrials; warning( ...
                    'Number of movies found  in the source folder is different from the number of trials.'); end

            % initialize
            % load batchprocess json (to have a list of init subhashes for the transformation outputs of each file)
            bias = 0; if contains(obj.folder,'2'); bias = numFiles; end
            filename = dir(fullfile(obj.folder,'*_init.json'));
            batchinit = readJson(fullfile(filename(1).folder,filename(1).name));
            initlist = batchinit.subHash(bias+(1:numFiles));

            % loop
            for i = 1:numFiles
                if sum(fail_idx(:,i))==0; continue; end % no fails for this trial: YAY!
                disp(fileList{i})

                % Load the movie data
                movie = robust_io('load',fileList2{i}, 'movie').movie;

                % load transformation maps
                filein = fullfile(obj.folder,'inits',[initlist{i},'.mat']);
                transform = robust_io('load',filein).process.operation;

                % list failed periods
                bps = convertPeriods(movie.badperiods(:,2:3),1);
                thisfails = convertPeriods(fail_idx(~bps,i)); % [start, end] without the bad periods
                numfails = height(thisfails);
                
                all_interp_tp = cell(numfails,1);
                for i_fail = 1:numfails
                    % extract transforms adjacent to the failed period
                    idx1 = thisfails(i_fail,1)-1;
                    idx2 = thisfails(i_fail,2)+1;
                    if idx1<1; idx1 = idx2; end
                    if idx2>size(transform.warp,4); idx2 = idx1; end
                    adj_tf = transform.warp(:,:,:,[idx1,idx2]);
                    nmissing = diff(thisfails(i_fail,:))+1;
                    interp_tp = zeros(movie.h,movie.w,2,nmissing+2);
                    
                    for i_dim = 1:2
                        cube = squeeze(adj_tf(:,:,i_dim,:));
                        [ny,nx,nz]=size(cube);
                        [x, y, z] = meshgrid(1:nx, 1:ny, 1:nz);
                        [xq, yq, zq] = meshgrid(1:nx, 1:ny, 1:1/(nmissing+1):nz);
                        interp_tp(:,:,i_dim,:) = interp3(x, y, z, cube, xq, yq, zq, 'spline');
                    end

                    transform.warp(:,:,:,thisfails(i_fail,1):thisfails(i_fail,2)) = interp_tp(:,:,:,2:end-1);

                    all_interp_tp{i_fail} = interp_tp;
                    
                end

                % build processor
                p = OpticFlowRegistration(movie);
                p.operation = transform;
                p = p.run;
                movie_corrected = BasicMovieProcessor('remove_badperiods',p.data_processed).run.data_processed;

                % Load the previously registred movie data
                movie_oldreg = robust_io('load',fileList{i}, 'movie').movie;
                movie_oldreg = BasicMovieProcessor('remove_badperiods',movie_oldreg).run.data_processed;

                if ~pl; continue; end

                for i_fail = 1:numfails
                    interp_tp = all_interp_tp{i_fail};

                    buffer = 20;


                    frcorrnew = avg_frame_correlation(movie_corrected.stack( ...
                        buffer:end-buffer,buffer:end-buffer, :));
                    frcorrold = avg_frame_correlation(movie_oldreg.stack( ...
                        buffer:end-buffer,buffer:end-buffer, :));
                    frcorrnew = nanzscore(frcorrnew);
                    frcorrold = nanzscore(frcorrold);

                    figure;
                    nrows = 3;
                    ncols = size(interp_tp,4)+1; % nmissing + adjacent frames + extra line plot
                    
                    for i_fr = 1:ncols-1
                        % get the actual frame numbers (but without
                        % counting the bad periods!) - so these are NOT
                        % necessarily the correct frame numbers in context
                        actualfrnum = thisfails(i_fail,1)+i_fr-2;

                        % real part
                        subplot(nrows,ncols,i_fr)
                        imagesc(interp_tp(:,:,1,i_fr)); axis square; xticklabels([]); yticklabels([]);
                        title(['fr #',num2str(actualfrnum)])
                        if i_fr == 1; ylabel('real'); end
    
                        % imaginary part
                        subplot(nrows,ncols,i_fr + ncols)
                        imagesc(interp_tp(:,:,2,i_fr)); axis square; xticklabels([]); yticklabels([]);
                        if i_fr == 1; ylabel('imaginary'); end
    
                        % corrected frame
                        subplot(nrows,ncols,i_fr + ncols*2)
                        imagesc(movie_corrected.stack(:,:,actualfrnum)); colormap('gray');
                        axis square; xticklabels([]); yticklabels([]);
                        if i_fr == 1; ylabel('adj. frames'); end
                    end

                    % restult correlation comparison
                    subplot(1,ncols,ncols)
                    overhang = 3; % including adjacent frames, one-sided
                    idx1 = max([1,thisfails(i_fail,1)-overhang]);
                    idx2 = min([movie.nfr,thisfails(i_fail,2)+overhang]);
                    y = [frcorrold((idx1:idx2)-1), frcorrnew((idx1:idx2)-1)];
                    b = plot(y,'LineWidth',2);
                    hold on
                    line([1 1]*overhang+.5 , ylim , ...
                        'Color','r', 'LineStyle', '--', 'LineWidth',1)
                    line([1 1]*height(y)-overhang+.5 , ylim , ...
                        'Color','r', 'LineStyle', '--', 'LineWidth',1)
                    legend(b, {'old warp','corrected'})
                    xticks(1:height(y))
                    xticklabels(idx1:idx2);
                    % ylim([-1 1])
                    ylabel('fr/fr z-scored correlation')
                    xlabel('frames')
                    axis tight
                end
            end

            cd(orig_folder)
        end

        function [fbf_corr, zsc_corr, putative_fails, hf] = QC(obj, th)
            % generate a frame-by-frame correlation curve

            orig_folder = fullfiletol(pwd);

            % get file list and move to subject-datapath
            [fileList, numFiles] = obj.sj.getFileListInSubfolder(obj.folder);
            if numFiles~=obj.sj.getNTrials; warning('Number of movies found is different from the number of trials.'); end

            % initialize
            fbf_corr = zeros(obj.sj.getNFrames-1, numFiles); % [numFrames-1, trials]

            % loop
            for i = 1:numFiles
                % Load the movie data
                disp(fileList{i})
                movie = robust_io('load',fileList{i}, 'movie').movie;

                buffer = 20; % registration algs can produce fast-variable 
                % edges (effect of warping at the edge of the image).
                % we take only the central portion of the image to estimate
                % the frame-by-frame correlation.
                fbf_corr(:,i) = avg_frame_correlation(movie.stack(buffer+1:end-buffer,buffer+1:end-buffer,:));
            end

            % get image size
            h = movie.h;
            w = movie.w;

            % zscored frame-by-frame correlation
            zsc_corr = nanzscore(fbf_corr);

            % define putatively failed frames
            putative_fails = [false(1,numFiles) ; zsc_corr < th]; % add one line at the start for frame 1
            % eliminate the last index of each putatively failed periods:
            % a good frame would still have low correlation with a bad
            % frame
            for i = 1:numFiles
                [thisfails, len] = convertPeriods(putative_fails(:,i));
                idx = diff(thisfails,[],2)>0; % everything should be positive in principle, but better have some robustness
                thisfails(idx,2) = thisfails(idx,2)-1;
                putative_fails(:,i) = convertPeriods(thisfails,true, len);
            end

            % return to original location
            cd(orig_folder)


            % PLOTTING

            hf = figure;
            hf.set('Position',[10 10 700 800]);


            ncol = 2;
            nrow = 3;
            [failsrows,failscols] = find(putative_fails);
            failscores = zsc_corr(putative_fails(2:end,:));
            [~, leastbadidx] = min(failscores,[],'omitmissing');
            [~, medianidx] = min(abs(failscores-median(failscores)),[],'omitmissing');
            [~, worstidx] = max(failscores,[],'omitmissing');

            subplot(nrow,ncol,1) 
            if ~isempty(leastbadidx)
                src_trial = failscols(leastbadidx);
                snip = Snippet(fileList{src_trial},failsrows(leastbadidx)-[1,0]);
                im = zeros(h,w,3);
                im(:,:,[1,3]) = snip.stack./max(snip.stack,[],'all','omitmissing');
                maxval = quantile(im(:),.97);
                imagesc(im./maxval); axis square
                xticks([]); yticks([])
                title(['least bad case: trial ',num2str(src_trial),', frame ',num2str(failsrows(leastbadidx))])
            end
            
            subplot(nrow,ncol,3) 
            if ~isempty(medianidx)
                src_trial = failscols(medianidx);
                snip = Snippet(fileList{src_trial},failsrows(medianidx)-[1,0]);
                im = zeros(h,w,3);
                im(:,:,[1,3]) = snip.stack./max(snip.stack,[],'all','omitmissing');
                maxval = quantile(im(:),.97);
                imagesc(im./maxval); axis square
                xticks([]); yticks([])
                title(['median case: trial ',num2str(src_trial),', frame ',num2str(failsrows(medianidx))])
            end

            subplot(nrow,ncol,5) 
            if ~isempty(worstidx)
                src_trial = failscols(worstidx);
                snip = Snippet(fileList{src_trial},failsrows(worstidx)-[1,0]);
                im = zeros(h,w,3);
                im(:,:,[1,3]) = snip.stack./max(snip.stack,[],'all','omitmissing');
                maxval = quantile(im(:),.97);
                imagesc(im./maxval); axis square
                xticks([]); yticks([])
                title(['worst case: trial ',num2str(src_trial),', frame ',num2str(failsrows(worstidx))])
            end

            subplot(nrow,ncol,2)
            imagesc(zsc_corr)
            title('z-scored frame correlations')
            xlabel('trial #')
            ylabel('frames')

            subplot(nrow,ncol,4)
            plot(zsc_corr); axis tight
            hold on
            bh = line(xlim,[th th],'Color','r','LineStyle','--');
            legend(bh,'thresh')
            title('z-scored frame correlations')
            xlabel('frames')


        end
    end
end