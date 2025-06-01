classdef PCAICA < MovieProcessing
    %%% comp on tiles works, but everything downstream assumes a single tile,
    %%% sorry!
    properties
        data_processed Movie
        operation
    end

    methods (Static)
        function movie_result = getFull2PMovieSubsampled(subject,folder)
            arguments
                subject Subject
                folder char
            end

            % knobs
            downsamplefactor = 2;
            
            % execute
            [fileList, numFiles] = subject.getFileListInSubfolder(folder);
            megabadperiods = [];
            for i = 1:numFiles
                disp(['Iter: #',num2str(i)])
                disp(['Loading ... ',fileList{i}])
                movie = robust_io('load',fileList{i}).movie; % loads 'movie' of type Movie
                
                % downsample in time
                movie = BasicMovieProcessor('downsamplet',movie).run('factor',downsamplefactor).data_processed;

                if i==1
                    totalnfr = movie.nfr*numFiles;
                    megastack = zeros(movie.h,movie.w,totalnfr);
                end
                idx = [1:movie.nfr] + movie.nfr*(i-1);
                thisbadperiods = movie.badperiods;
                thisbadperiods(:,2:3) = thisbadperiods(:,2:3) + movie.nfr*(i-1);
                megastack(:,:,idx) = movie.stack;
                megabadperiods = [megabadperiods; thisbadperiods];

                disp('//////////////')
            end
            movie_result = movie;
            movie_result.badperiods = megabadperiods;
            movie_result.stack = megastack;
            disp('done.')
        end

        function lin_crop = addnanpxsbackin(lin_crop, idx_nanpx)
            if isempty(idx_nanpx)
                return
            end
        
            tmp = nan(height(lin_crop)+numel(idx_nanpx),width(lin_crop));
            idx_numpx = true(height(lin_crop)+numel(idx_nanpx),1);
            idx_numpx(idx_nanpx) = false;
            tmp(idx_numpx,:) = lin_crop;
            lin_crop = tmp;
            clear tmp
        end

        function [lin_movie, idx_nanfr, idx_nanpx] = linearizeFrames(movie, elim_nanpx, idx_nanfr, idx_nanpx, method)
            arguments
                movie double
                elim_nanpx logical = false
                idx_nanfr = []
                idx_nanpx = []
                method char = 'svds'
            end

            % linearize frames (spatial structure is not relevant)
            [h,w,nfr] = size(movie);
            lin_movie = zeros(h*w,nfr);
            for i = 1:nfr
                thisframe = movie(:,:,i);
                lin_movie(:,i) = thisframe(:);
            end

            % temporarily take out any cols (frames) which are all nan
            if isempty(idx_nanfr)
                idx_nanfr = find(sum(isnan(lin_movie),1)==height(lin_movie));
                if ~isempty(idx_nanfr)
                    fprintf('Found following unexpected NaN frames: %s\n',num2str(idx_nanfr))
                end
            end
            lin_movie(:,idx_nanfr) = [];
            
            % eliminate any rows (px) with nan values (but keep the ids)
            if isempty(idx_nanpx)
                idx_nanpx = find(sum(isnan(lin_movie),2)>0);
                str = 'not ';
                if elim_nanpx; lin_movie(idx_nanpx,:) = []; str = ''; end
                if ~isempty(idx_nanpx)
                    fprintf('Found %s pixels containing NaNs: %seliminated.\n', ...
                        num2str(numel(idx_nanpx)),str)
                end
            elseif elim_nanpx
                lin_movie(idx_nanpx,:) = [];
            end
            
            % add nan frames back in
            if isempty(idx_nanfr) || strcmp(method,'svds'); return; end
            tmp = nan(height(lin_movie),width(lin_movie)+numel(idx_nanfr));
            idx_numfr = true(width(lin_movie)+numel(idx_nanfr),1);
            idx_numfr(idx_nanfr) = false;
            tmp(:,idx_numfr) = lin_movie;
            lin_movie = tmp;
            clear tmp
        end
    end

    methods %(Access=protected)
        function [ICAresults, params, hf] = spatialICA(obj, PCAresults, numcomponentsICA)
            arguments
                obj
                PCAresults struct
                numcomponentsICA double = 10; % [min: 2, max: numcomponentsPCA]
            end
            
            % init
            numcomponentsPCA = PCAresults.numcomponents;
            score = PCAresults.score;
            
            fprintf('Running ICA on the first %s PCs to obtain %s ICs ...\n',...
                num2str(numcomponentsPCA), num2str(numcomponentsICA)); tic
            
            Mdl = rica(score,numcomponentsICA);
            
            data_ica = transform(Mdl, score);
            
            % plot Independent components
            hf = figure;
            subplot(131); plot(data_ica.^2); ylabel('squared ic measurements'); xlabel('frames')
            subplot(132); semilogy(data_ica.^2); xlabel('frames')
            subplot(133); imagesc(data_ica); xlabel('ICs'); ylabel('frames')
            
            toc
            
            % return
            ICAresults.data = data_ica;
            ICAresults.Mdl = Mdl;
            
            params.numcomponentsPCA = numcomponentsPCA;
            params.numcomponentsICA = numcomponentsICA;
        end

        function [PCAresults, reconstructions, params, hf] = spatialPCA(obj,data_raw)
            arguments
                obj
                data_raw Movie
            end
            
            % init
            crop = data_raw.stack;
            nfr = data_raw.nfr;
            
            % init default parameters
            params = [];
            
            % linearize movie frames and eliminate NaNs
            [lin_crop, idx_nanfr, idx_nanpx] = obj.linearizeFrames(crop, true,[],[],obj.init.method);
            
            %% define dataset
            
            a = lin_crop'; % [time, px]
            [a,mu,stda] = nanzscore(a);
            
            %% cropped pca
            disp('Running PCA ...'); tic
            
            if strcmp(obj.init.method,'svds')
                k = 50;  % number of components to compute
                [U, S, V] = svds(a, k);  % a' is (pixels x time)
                score = U * S;       % projections (like PCA scores)
                coef = V;                % principal directions (time)
                latent_all = diag(S).^2; % eigenvalues
                explained = latent_all / sum(latent_all);
            else
                % [coef,score,~,~,explained,~] = pca(a,'Economy',false,'Algorithm','als'); % vars are px
                [coef,~,~,~,explained,~] = pca(a,'Economy',false,'Rows','pairwise'); % vars are px
            end
            numcomponentsPCA = find(cumsum(explained)>.99,1); % explain 99% of variance
            fullscore = a*coef; score = fullscore(:,1:numcomponentsPCA);

            figure; subplot(121); imagesc(fullscore), title('score'); subplot(122); imagesc(coef), title('coef');
            figure; plot(explained); xlabel('PC #'); ylabel('explained variance')
            figure; plot(cumsum(explained)); xlabel('PC #'); ylabel('cum. explained variance')

            toc
            
            %% (optional) reconstruct from PCA
            
            lin_crop_rc = transpose(fullscore*coef');
            
            % lin_crop_fullyrc = lin_crop_rc.*stda' + repmat(mu',1,nfr);
            lin_crop_fullyrc = lin_crop_rc.*stda' + mu';
            
            % add pixels containing nans back in (this time as nans only)
            % add nan frames back in
            lin_crop = obj.addnanpxsbackin(lin_crop, idx_nanpx);
            lin_crop_rc = obj.addnanpxsbackin(lin_crop_rc, idx_nanpx);
            lin_crop_fullyrc = obj.addnanpxsbackin(lin_crop_fullyrc, idx_nanpx);
            
            % reconstruct 3d movie
            crop_rc = reshape(lin_crop_rc, size(crop));
            crop_fullyrc = reshape(lin_crop_fullyrc, size(crop));
            
            % plot
            hf = figure;
            subplot(231); imagesc(Movie(crop).timeavg); title('raw')
            subplot(232); imagesc(Movie(crop_rc).timeavg); title('reconstructed')
            subplot(233); imagesc(Movie(crop_fullyrc).timeavg); title('fully reconstr')
            subplot(234); imagesc(lin_crop)
            subplot(235); imagesc(lin_crop_rc)
            subplot(236); imagesc(lin_crop_fullyrc)
            
            %% store to structures
            
            % results
            PCAresults.rawinput = crop; % input data
            PCAresults.linearizedinput = lin_crop;
            PCAresults.zscoredinput = a;
            PCAresults.inputstd = stda;
            PCAresults.inputmu = mu;
            PCAresults.fullscore = fullscore;
            PCAresults.score = score;
            PCAresults.coef = coef;
            PCAresults.numcomponents = numcomponentsPCA;
            PCAresults.idx_nanpx = idx_nanpx;
            PCAresults.idx_nanfr = idx_nanfr; % these are the idx of 
            % nan-containing frames that were taken out AFTER removing 
            % recorded 'badperiods'
            
            % reconstructed stacks
            reconstructions.raw = crop; % input data
            reconstructions.zscore = crop_rc;
            reconstructions.full = crop_fullyrc;
            
            % parameters in params
        end

        function results = run_singletile(obj, stack)            
            crop = Movie(stack);
            figure; imagesc(crop.timeavg); colormap('gray')
            
            % PCA
            [PCAresults,~,~,~] = obj.spatialPCA(crop);
            
            % ICA
            [ICAresults,~,~] = obj.spatialICA(PCAresults, PCAresults.numcomponents);

            % return
            results.PCA = PCAresults;
            results.ICA = ICAresults;

        end

        function [results, tile_coords] = apply_to_tiles(obj, stack, myFunction)
        % Applies myFunction to each spatial tile in stack (Y, X, T)
        % Returns a cell array of outputs and tile coordinates

            tile_size = obj.operation.tile_size;
            tile_overlap = obj.operation.tile_overlap;
        
            [H, W, ~] = size(stack);
            tx = tile_size(1);
            ty = tile_size(2);
            ox = round(tile_overlap(1) * tx);
            oy = round(tile_overlap(2) * ty);
        
            % Define tile positions
            x_starts = 1 : (tx - ox) : max(W - tx + 1, 1);
            y_starts = 1 : (ty - oy) : max(H - ty + 1, 1);
        
            n_tiles = length(x_starts) * length(y_starts);
            results = cell(n_tiles, 1);
            tile_coords = cell(n_tiles, 1);
        
            tile_idx = 1;
            for yi = 1:length(y_starts)
                for xi = 1:length(x_starts)
                    x0 = x_starts(xi);
                    y0 = y_starts(yi);
                    x1 = min(x0 + tx - 1, W);
                    y1 = min(y0 + ty - 1, H);
        
                    tile = stack(y0:y1, x0:x1, :);
        
                    % Apply user function
                    results{tile_idx} = myFunction(tile);
                    tile_coords{tile_idx} = [x0, y0; x1, y1];
        
                    tile_idx = tile_idx + 1;
                end
            end
        end

    end

    methods
        function obj = PCAICA(data_raw)
            arguments
                data_raw Movie = Movie('')
            end

            obj = obj.setDataRaw(data_raw);
            
            % set default operation params
            obj.operation.tile_size = [600 600]; % [px] (for actual tiling try [170 170])
            obj.operation.tile_overlap = [.5 .5]; % fraction of tile overlap
            obj.operation.th_val = 10;
            obj.operation.th_occupancy = .5; % set =nan only if the px is <th_val at least 50% of the time
            obj.operation.blursigma = 0; % use value of 1 to actually blur a little but not too much
            obj.operation.rollavg_win = 0; % use value of 3 to actually blur a little but not too much
            obj.operation.downsamplefactor = 2;
            obj.operation.resize_factor = 2;

            % set default method
            obj.init.method = 'svds';
        end

        function [obj, results, tile_coords, pcaica_w, ROImap, hf, stack] = run(obj, ROImap)
            % header
            obj.disp_runheader;

            % preprocessing steps
            stack = obj.preprocessStack;
            
            % run PCAICA on tiles
            [results, tile_coords] = obj.apply_to_tiles(stack, @obj.run_singletile);
            obj.operation.results = results;
            obj.operation.tile_coords = tile_coords;

            % save results
            obj.data_processed = Movie('');
            s.pcaica = obj.getLight;
            robust_io('save','pcaica.mat',s,'-v7.3'); clear s
            
            % impact of individual ICs (assumes a single tile)
            if ~exist("ROImap",'var'); ROImap = []; end
            [ROImap, hf, spectral, wmaps] = plotICimpact(obj, results{1}.ICA, results{1}.PCA, [], ROImap, false);

            % ultralight save (contains some spectral analysis results)
            pcaica_w.PCcoef = results{1}.PCA.coef;
            pcaica_w.ICs = results{1}.ICA.data;
            pcaica_w.idx_nanpx = results{1}.PCA.idx_nanpx;
            pcaica_w.idx_nanfr = results{1}.PCA.idx_nanfr;
            pcaica_w.Mdl = results{1}.ICA.Mdl;
            pcaica_w.ICweightmaps = wmaps;
            pcaica_w.pxstd = results{1}.PCA.inputstd;
            [~,pcaica_w.pxmufullsize,pcaica_w.pxstdfullsize] = ...
                nanzscore(obj.linearizeFrames(obj.data_raw.stack,true,[],[],obj.init.method)');
            pcaica_w.pxmu = results{1}.PCA.inputmu;
            pcaica_w.numusedPCs = results{1}.PCA.numcomponents;
            pcaica_w.spectral_bias = spectral;

            s.pcaica_w = pcaica_w;
            robust_io('save','pcaica_w.mat',s,'-v7.3');
            
            save_figures_to_pdf(hf,mfilename);
        end

        function stack = preprocessStack(obj)
            movie = obj.data_raw;

            % extract free params
            th_val = obj.operation.th_val;
            th_occupancy = obj.operation.th_occupancy;
            blursigma = obj.operation.blursigma;
            rollavg_win = obj.operation.rollavg_win;
            downsamplefactor = obj.operation.downsamplefactor;
            resize_factor = obj.operation.resize_factor;

            % prepare movie
            
            % remove bad periods
            movie_clean = BasicMovieProcessor('remove_badperiods',movie).run.data_processed;
            
            % rolling average in time
            movie_clean = BasicMovieProcessor('movmean',movie_clean).run('win_size',rollavg_win).data_processed;

            % downsample in time (commented out because done in getFull2PMovieSubsampled() already)
            % movie_clean = BasicMovieProcessor('downsamplet',movie_clean).run('factor',downsamplefactor).data_processed;
            
            
            % init
            stack = movie_clean.stack;
            h = movie_clean.h;
            w = movie_clean.w;
            nfr = movie_clean.nfr;
            
            % set all values under a threshold = nan (to be done before any spatial operations!)
            [rows,cols] = find(sum(stack<th_val,3)./nfr >= th_occupancy);
            for i = 1:length(rows)
                stack(rows(i), cols(i), :) = NaN;
            end
            
            % apply 2d gaussian blur
            stack = BasicMovieProcessor('gauss_blur2d',Movie(stack)).run('sigma',blursigma).data_processed.stack;

            % resize image (downsample based on nearest-neighbor interpolation)
            % Method 'nearest' is chosen because it leads to the highest
            % resconstruction accuracy when computing 
            % imresize(scale) -> imresize(1/scale)
            stack = imresize(stack,'nearest','scale',1/resize_factor);
        end

        function [ROImap, hf, spectral_bias, wmapsout] = plotICimpact(obj, ICAresults, PCAresults, idx,ROImap, keepdrawing) %% TENTATIVELY FIXED
            arguments
                obj
                ICAresults struct
                PCAresults struct
                idx double = []
                ROImap double = []
                keepdrawing logical = false
            end

            % init vars
            data_ica = ICAresults.data;
            Mdl = ICAresults.Mdl;
            numcomponentsICA = width(data_ica);
            mu = PCAresults.inputmu;
            stda = PCAresults.inputstd;
            coef = PCAresults.coef;
            numcomponentsPCA = PCAresults.numcomponents;
            crop = PCAresults.rawinput;
            lin_crop = PCAresults.linearizedinput;
            nfr = size(crop,3); % equivalent to Movie(crop).nfr but faster
            idx_nanpx = PCAresults.idx_nanpx;
            
            if ~exist('idx','var') || isempty(idx)
                [~, idx] = max(var(data_ica)); % pick most variant IC
            end

            %% analyse spectral biases of IC weight maps

            % get weight maps
            wmaps = cell(numcomponentsICA,1);
            for i = 1:numcomponentsICA
                weightmap = Mdl.TransformWeights(:,i)'*coef(:,1:numcomponentsPCA)';
                weightmap = obj.addnanpxsbackin(weightmap', idx_nanpx);
                weightmap = reshape(weightmap,height(crop),width(crop));

                wmaps{i} = weightmap;
            end
            wmapsout.whole = wmaps;


            % crop weight maps (excude borders, correlation drivers because of alignment artefacts)
            buffer = 20;
            wmaps = cellfun(@(x) x(buffer+1:end-buffer,buffer+1:end-buffer),wmaps,'UniformOutput',false);
            wmapsout.cropped = wmaps;

            % index by decreasing IC variance in time
            [~, varidx] = sort(var(data_ica), 'descend');
            
            hf = gobjects(1,10);

            % plotting weight maps
            hf(1) = figure;
            hf(2) = figure;
            ncols = floor(numcomponentsICA/2)+1;
            nrows = floor(numcomponentsICA/(ncols-1));
            t = tiledlayout(nrows, ncols, 'TileSpacing', 'compact', 'Padding', 'compact');
            for i = 1:numcomponentsICA
                nexttile
                weightmap = wmaps{varidx(i)}; % plot in order of decreasing IC variance in time

                figure(hf(1));
                %subplot(2,5,i)
                y = weightmap;
                imagesc(y)
                title(['IC #',num2str(varidx(i)),', var=',num2str(var(data_ica(:,varidx(i))),'%.1f')])
                axis square
                xticklabels([]),yticklabels([])

                figure(hf(2));
                %subplot(2,5,i)
                y = BasicMovieProcessor('clahe',Movie(weightmap)).run().data_processed.timeavg;
                imagesc(y)
                title(['IC #',num2str(varidx(i)),', var=',num2str(var(data_ica(:,varidx(i))),'%.1f')])
                axis square
                xticklabels([]),yticklabels([])
            end

            % extract spectral bias metrics for IC weight maps
            [results, hf_tmp] = sortImgsPOW(wmaps, true);
            for i_fig = 1:numel(hf_tmp)
                hf(i_fig+2) = hf_tmp(i_fig);
            end
            clear hf_tmp

            % plot IC variance in time vs IC weight map spectral bias score
            hf(7) = figure;
            x = std(data_ica);
            y = results.spectral_bias.bias_scores;
            scatter(x,y, 'filled'); axis square
            xlabel('IC STD')
            ylabel('Spectral Bias Score')
            hold on
            rsquared = fitlm(x, y).Rsquared.Ordinary;
            text(min(xlim)+.2*diff(xlim),mean(y), ...
                ['R^2 = ',num2str(rsquared)],'FontSize',8,'Color','r');
            for i = 1:numcomponentsICA % Add index labels
                text(x(i) + 0.03*diff(xlim), y(i), num2str(i), 'FontSize', 8,'Color','b');
            end
            set(gcf,"Position",[1,1,300,200])

            hf(8) = figure;
            [~,biasidx] = sort(results.spectral_bias.bias_scores,'descend');
            for i = 1:numcomponentsICA
                nexttile
                weightmap = wmaps{biasidx(i)}; % plot in order of decreasing IC variance in time

                %subplot(2,5,i)
                y = BasicMovieProcessor('clahe',Movie(weightmap)).run().data_processed.timeavg;
                imagesc(y)
                title(['IC #',num2str(biasidx(i)),', var=',num2str(var(data_ica(:,biasidx(i))),'%.1f')])
                axis square
                xticklabels([]),yticklabels([])
            end

            spectral_bias = results.spectral_bias;
            spectral_bias.ICtimevar = var(data_ica)';


            %% zero out ICs
            
            % keep only idx
            to_elim = true(width(data_ica),1); to_elim(idx) = false;
            
            data_ica_clean = data_ica;
            data_ica_clean(:,to_elim) = zeros(height(data_ica),sum(to_elim));
            
            %% reconstruct data
            
            % inverse ICA
            data_pca = data_ica_clean*Mdl.TransformWeights';
            
            % inverce PCA
            lin_crop_rc = transpose(data_pca*coef(:,1:numcomponentsPCA)');
            lin_crop_fullyrc = lin_crop_rc.*stda';% + mu';
            
            % add pixels containing nans back in (this time as nans only)
            % add nan frames back in
            lin_crop = obj.addnanpxsbackin(lin_crop, idx_nanpx);
            lin_crop_rc = obj.addnanpxsbackin(lin_crop_rc, idx_nanpx);
            lin_crop_fullyrc = obj.addnanpxsbackin(lin_crop_fullyrc, idx_nanpx);
            
            % reconstruct 3d movie
            crop_rc = reshape(lin_crop_rc, size(crop));
            crop_fullyrc = reshape(lin_crop_fullyrc, size(crop));

            crop_rc = crop - crop_rc;
            crop_fullyrc = crop - crop_fullyrc;
            
            hf(3) = figure;
            subplot(231); imagesc(Movie(crop).timeavg); title('raw')
            subplot(232); imagesc(Movie(crop_rc).timeavg); title('reconstructed')
            subplot(233); imagesc(Movie(crop_fullyrc).timeavg); title('fully reconstr')
            subplot(234); imagesc(lin_crop)
            subplot(235); imagesc(lin_crop_rc)
            subplot(236); imagesc(lin_crop_fullyrc)
            
            %% activity of ROIs

            return; % ############
            
            if isempty(ROImap) || keepdrawing
                ROImap = imageSequenceGUI({Movie(crop).timeavg},{computeLocalCorrelationMap(crop)},ROImap,'select some ROIs');
            end
            rois = unique(ROImap);
            
            activity = zeros(nfr, numel(rois));
            for i_fr = 1:nfr
                thisfr = crop(:,:,i_fr);
                for i_roi = 1:numel(rois)
                    activity(i_fr,i_roi) = mean(thisfr(ROImap==rois(i_roi)),'all','omitmissing');
                end
            end
            activity_rc = zeros(nfr, numel(rois));
            for i_fr = 1:nfr
                thisfr = crop_rc(:,:,i_fr);
                for i_roi = 1:numel(rois)
                    activity_rc(i_fr,i_roi) = mean(thisfr(ROImap==rois(i_roi)),'all','omitmissing');
                end
            end
            activity_fullyrc = zeros(nfr, numel(rois));
            for i_fr = 1:nfr
                thisfr = crop_fullyrc(:,:,i_fr);
                for i_roi = 1:numel(rois)
                    activity_fullyrc(i_fr,i_roi) = mean(thisfr(ROImap==rois(i_roi)),'all','omitmissing');
                end
            end
            
            hf(9) = figure;
            subplot(141); plot(activity); xlabel('frames'); ylabel('ROI activity (raw)'); ax(1) = gca;
            labels = [{'background'}, arrayfun(@(i) sprintf('cell %d', i), 1:(numel(rois)-2), 'UniformOutput', false), {'cell-sized bg'}];
            legend(labels)
            subplot(142); plot(activity_rc); xlabel('frames'); ylabel('ROI activity (reconst)'); ax(2) = gca;
            subplot(143); plot(activity_fullyrc); xlabel('frames'); ylabel('ROI activity (fully reconst)'); ax(3) = gca;
            subplot(144); plot(activity - repmat(activity(:,1),1,width(activity)));ax(4) = gca;
            xlabel('frames'); ylabel('ROI activity (raw - background)')
            linkaxes(ax,'x')

            hf(10) = figure;
            subplot(141); imagesc(activity); ax(1)=gca;
            subplot(142); imagesc(activity_rc); ax(2)=gca;
            subplot(143); imagesc(activity_fullyrc); ax(3)=gca;
            subplot(144); imagesc(activity - repmat(activity(:,1),1,width(activity))); ax(4)=gca;
            linkaxes(ax,'xy')
        end


    end
end