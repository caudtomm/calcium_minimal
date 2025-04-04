classdef PCAICA < MovieProcessing
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
                    load(fileList{i}) % loads 'movie' of type Movie
                    
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
    end

    methods (Static, Access=protected)
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
    end

    methods (Access=protected)
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
            
            %% linearize frames (spatial structure is not relevant)
            lin_crop = zeros(height(crop)*width(crop),nfr);
            for i = 1:nfr
                thisframe = crop(:,:,i);
                lin_crop(:,i) = thisframe(:);
            end
            
            % temporarily take out any cols (frames) which are all nan
            idx_nanfr = find(sum(isnan(lin_crop),1)==height(lin_crop));
            if ~isempty(idx_nanfr)
                fprintf('Found following unexpected NaN frames: %s\n',num2str(idx_nanfr))
            end
            lin_crop(:,idx_nanfr) = [];
            
            % eliminate any rows (px) with nan values (but keep the ids)
            idx_nanpx = find(sum(isnan(lin_crop),2)>0);
            lin_crop(idx_nanpx,:) = [];
            if ~isempty(idx_nanpx)
                fprintf('Found %s pixels containing NaNs: eliminated.\n',num2str(numel(idx_nanpx)))
            end
            
            % add nan frames back in
            if ~isempty(idx_nanfr)
                tmp = nan(height(lin_crop),width(lin_crop)+numel(idx_nanfr));
                idx_numfr = true(width(lin_crop)+numel(idx_nanfr),1);
                idx_numfr(idx_nanfr) = false;
                tmp(:,idx_numfr) = lin_crop;
                lin_crop = tmp;
                clear tmp
            end
            
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
            numcomponentsPCA = find(cumsum(explained)>.95,1); % explain 95% of variance
            fullscore = a*coef; score = fullscore(:,1:numcomponentsPCA);

            figure; subplot(121); imagesc(fullscore), title('score'); subplot(122); imagesc(coef), title('coef');
            
            toc
            
            %% (optional) reconstruct from PCA
            
            lin_crop_rc = transpose(fullscore*coef');
            
            lin_crop_fullyrc = lin_crop_rc.*stda' + repmat(mu',1,nfr);
            
            % add pixels containing nans back in (this time as nans only)
            % add nan frames back in
            lin_crop = obj.addnanpxsbackin(lin_crop, idx_nanpx);
            lin_crop_rc = obj.addnanpxsbackin(lin_crop_rc, idx_nanpx);
            lin_crop_fullyrc = obj.addnanpxsbackin(lin_crop_fullyrc, idx_nanpx);
            
            crop_rc = zeros(size(crop));
            for i = 1:nfr
                crop_rc(:,:,i) = reshape(lin_crop_rc(:,i),[height(crop),width(crop)]);
            end
            
            crop_fullyrc = zeros(size(crop));
            for i = 1:nfr
                crop_fullyrc(:,:,i) = reshape(lin_crop_fullyrc(:,i),[height(crop),width(crop)]);
            end
            
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
            [ICAresults,~,~] = obj.spatialICA(PCAresults, 10);

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
            obj.operation.tile_size = [170 170]; % [px]
            obj.operation.tile_overlap = [.5 .5]; % fraction of tile overlap
            obj.operation.th_val = 10;
            obj.operation.th_occupancy = .5; % set =nan only if the px is <th_val at least 50% of the time
            obj.operation.blursigma = 1;
            obj.operation.rollavg_win = 3;
            obj.operation.downsamplefactor = 2;
        end

        function [obj, results, tile_coords] = run(obj)
            % header
            obj.disp_runheader;

            % preprocessing steps
            stack = obj.preprocessStack;
            
            % run PCAICA on tiles
            [results, tile_coords] = obj.apply_to_tiles(stack, @obj.run_singletile);
            
            % impact of individual ICs
            

            % return

            % tail
            % obj = obj.tailsequence;
        end

        function stack = preprocessStack(obj)
            movie = obj.data_raw;

            % extract free params
            tile_size = obj.operation.tile_size;
            tile_overlap = obj.operation.tile_overlap;
            th_val = obj.operation.th_val;
            th_occupancy = obj.operation.th_occupancy;
            blursigma = obj.operation.blursigma;
            rollavg_win = obj.operation.rollavg_win;
            downsamplefactor = obj.operation.downsamplefactor;

            % prepare movie
            
            % remove bad periods
            movie_clean = BasicMovieProcessor('remove_badperiods',movie).run.data_processed;
            
            % rolling average in time
            % movie_clean = BasicMovieProcessor('movmean',movie_clean).run('win_size',rollavg_win).data_processed;
            % 
            % % downsample in time
            % movie_clean = BasicMovieProcessor('downsamplet',movie_clean).run('factor',downsamplefactor).data_processed;
            
            
            % init
            stack = movie_clean.stack;
            h = movie_clean.h;
            w = movie_clean.w;
            nfr = movie_clean.nfr;
            
            % set all values under a threshold = nan
            [rows,cols] = find(sum(stack<th_val,3)./nfr >= th_occupancy);
            for i = 1:length(rows)
                stack(rows(i), cols(i), :) = NaN;
            end
            
            % apply 2d gaussian blur
            % stack = BasicMovieProcessor('gauss_blur2d',Movie(stack)).run('sigma',blursigma).data_processed.stack;
        end


        function [ROImap, hf] = plotICimpact(obj, ICAresults, PCAresults, idx) %% TENTATIVELY FIXED
            % init vars
            data_ica = ICAresults.data;
            Mdl = ICAresults.Mdl;
            mu = PCAresults.inputmu;
            stda = PCAresults.inputstd;
            coef = PCAresults.coef;
            numcomponentsPCA = PCAresults.numcomponents;
            crop = PCAresults.rawinput;
            lin_crop = PCAresults.linearizedinput;
            nfr = obj.data_raw.nfr;
            
            if ~exist('idx','var') || isempty(idx)
                [~, idx] = max(var(data_ica)); % pick most variant IC
            end
            
            %% zero out ICs
            
            % keep only idx
            to_elim = true(width(data_ica),1); to_elim(idx) = false;
            
            % elim only idx
            to_elim = false(width(data_ica),1); to_elim(idx) = true;
            
            data_ica_clean = data_ica;
            data_ica_clean(:,to_elim) = zeros(height(data_ica),sum(to_elim));
            
            %% reconstruct data
            
            % inverse ICA
            data_pca = data_ica_clean*Mdl.TransformWeights';
            
            % inverce PCA
            lin_crop_rc = transpose(data_pca*coef(:,1:numcomponentsPCA)');
            lin_crop_fullyrc = lin_crop_rc.*stda' + repmat(mu',1,nfr);
            
            % add pixels containing nans back in (this time as nans only)
            % add nan frames back in
            if ~isempty(idx_nanpx)
                % lin_crop
                tmp = nan(height(lin_crop)+numel(idx_nanpx),width(lin_crop));
                idx_numpx = true(height(lin_crop)+numel(idx_nanpx),1);
                idx_numpx(idx_nanpx) = false;
                tmp(idx_numpx,:) = lin_crop;
                lin_crop = tmp;
                clear tmp
                
                % lin_crop_rc
                tmp = nan(height(lin_crop_rc)+numel(idx_nanpx),width(lin_crop_rc));
                idx_numpx = true(height(lin_crop_rc)+numel(idx_nanpx),1);
                idx_numpx(idx_nanpx) = false;
                tmp(idx_numpx,:) = lin_crop_rc;
                lin_crop_rc = tmp;
                clear tmp
                
                % lin_crop_fullyrc
                tmp = nan(height(lin_crop_fullyrc)+numel(idx_nanpx),width(lin_crop_fullyrc));
                idx_numpx = true(height(lin_crop_fullyrc)+numel(idx_nanpx),1);
                idx_numpx(idx_nanpx) = false;
                tmp(idx_numpx,:) = lin_crop_fullyrc;
                lin_crop_fullyrc = tmp;
                clear tmp
            end
            
            % reconstruct 3d movie
            crop_rc = zeros(size(crop));
            for i = 1:nfr
                crop_rc(:,:,i) = reshape(lin_crop_rc(:,i),[height(crop),width(crop)]);
            end
            
            crop_fullyrc = zeros(size(crop));
            for i = 1:nfr
                crop_fullyrc(:,:,i) = reshape(lin_crop_fullyrc(:,i),[height(crop),width(crop)]);
            end
            
            hf(1) = figure;
            subplot(231); imagesc(Movie(crop).timeavg); title('raw')
            subplot(232); imagesc(Movie(crop_rc).timeavg); title('reconstructed')
            subplot(233); imagesc(Movie(crop_fullyrc).timeavg); title('fully reconstr')
            subplot(234); imagesc(lin_crop)
            subplot(235); imagesc(lin_crop_rc)
            subplot(236); imagesc(lin_crop_fullyrc)
            
            %% activity of ROIs
             
            ROImap = imageSequenceGUI({Movie(crop).timeavg},{computeLocalCorrelationMap(crop)},ROImap,'select some ROIs');
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
            
            hf(2) = figure;
            subplot(141); plot(activity); xlabel('frames'); ylabel('ROI activity (raw)');
            labels = [{'background'}, arrayfun(@(i) sprintf('cell %d', i), 1:(numel(rois)-2), 'UniformOutput', false), {'cell-sized bg'}];
            legend(labels)
            subplot(142); plot(activity_rc); xlabel('frames'); ylabel('ROI activity (reconst)')
            subplot(143); plot(activity_fullyrc); xlabel('frames'); ylabel('ROI activity (fully reconst)')
            subplot(144); plot(activity - repmat(activity(:,1),1,width(activity)));
            xlabel('frames'); ylabel('ROI activity (raw - background)')
            
        end
    end
end