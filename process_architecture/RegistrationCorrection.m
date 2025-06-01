classdef RegistrationCorrection < MovieProcessing
    properties
        data_processed
        operation
    end

    methods (Static)
        function [hf, frcorrnew, frcorrold] = auditSingleCorrection(thisfail, interp_tp, movie_corrected, movie_oldreg, pl)
            % interp_tp = all_interp_tp{i_fail};
            % thisfail = thisfails(i_fail,:);

            hf = [];

            buffer = 20;

            frcorrnew = avg_frame_correlation(movie_corrected.stack( ...
                buffer:end-buffer,buffer:end-buffer, :));
            frcorrold = avg_frame_correlation(movie_oldreg.stack( ...
                buffer:end-buffer,buffer:end-buffer, :));
            frcorrnew = nanzscore(frcorrnew);
            frcorrold = nanzscore(frcorrold);

            if ~pl; return; end

            hf = figure;
            nrows = 3;
            ncols = size(interp_tp,4)+1; % nmissing + adjacent frames + extra line plot
            
            for i_fr = 1:ncols-1
                % get the actual frame numbers (but without
                % counting the bad periods!) - so these are NOT
                % necessarily the correct frame numbers in context
                actualfrnum = thisfail(1)+i_fr-2;

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
            idx1 = max([1,thisfail(1)-overhang]);
            idx2 = min([movie_oldreg.nfr,thisfail(2)+overhang]);
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

    methods
        function obj = RegistrationCorrection(movie_in)
            arguments
                movie_in = ''
            end
            obj = obj@MovieProcessing(movie_in);

            % initialization
            obj.data_processed = obj.data_raw;
        end

        function obj = run(obj)
            % init
            fail_idx = obj.operation.fail_idx; % column vector
            transforminit = obj.operation.transforminit; % char array
            th = obj.operation.zscorethresh;
            movie = obj.data_raw; % Movie
            preregistered_movie = obj.operation.preregmovie;
            obj.data_processed = obj.data_raw;

            if sum(fail_idx)==0 % 
                disp('It looks like there are no fails for this trial: YAY!');
                movie_oldreg = robust_io('load',preregistered_movie, 'movie').movie;
                obj.data_processed = movie_oldreg;
                return
            end

            % load transformation maps;
            transform = robust_io('load',transforminit).process.operation;

            % list failed periods
            bps = convertPeriods(movie.badperiods(:,2:3),1,movie.nfr);
            thisfails = convertPeriods(fail_idx(~bps)); % [start, end] without the bad periods
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
            % movie_oldreg = robust_io('load',preregistered_movie, 'movie').movie;
            % movie_oldreg = BasicMovieProcessor('remove_badperiods',movie_oldreg).run.data_processed;

            % mark correction as good or bad
            buffer = 20; % spatial buffer (avoid contamination from registration borders)
            frcorrnew = avg_frame_correlation(movie_corrected.stack( ...
            buffer:end-buffer,buffer:end-buffer, :)); 
            frcorrnew = [0; nanzscore(frcorrnew)]; % nfr preserved this way
            % put into full lenght double array with nans where orig badperiods
            tmp = nan(movie.nfr,1); tmp(~bps) = frcorrnew; % this is a bit brittle. we expect sum(~bp) = length(frcorrnew), but this is not necessarily true if
            frcorrnew = tmp; clear tmp
            idx_bad = frcorrnew<th;
            bps(idx_bad) = true;
            bps = convertPeriods(bps);

            firstcol = [];
            if ~isempty(movie.badperiods); firstcol = movie.badperiods(1); end
            if isempty(firstcol); try; firstcol = movie.path.trial_num; catch; end; end
            if isempty(firstcol); firstcol = 1; end
            bps = [firstcol*ones(height(bps),1) , bps];
            
            movie_corrected = p.data_processed;
            movie_corrected.badperiods = bps;

            obj.data_processed = BasicMovieProcessor('replace_badperiods_with_nans',movie_corrected).run.data_processed;

            
        end
    end

end