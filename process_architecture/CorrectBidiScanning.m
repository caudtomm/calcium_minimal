classdef CorrectBidiScanning < MovieProcessing
    % - compute correction based on reference image 
    % (in which case save figure to a png and matlab fig) 
    % or load correction factor
    % - make sure correction factor is saved
    % - run applies precomputed correction to input movie
properties
    data_processed Movie    % inherited from MovieProcessing
    operation               % inherited from SingleProcess
end

properties (SetAccess=protected)
    reference_img double
end

methods (Static)
    function [images, phase_shifts_all, ...
            g_consensus] = ...
            estimate_g_and_correct_multiple(S, num_iters)
        % S is a 3D matrix of size (N, num_pairs, 2), where
        %   - S(:, i, 1) = s1 for pair i
        %   - S(:, i, 2) = s2 for pair i
        % num_iters: number of refinement iterations
        
        [N, num_pairs, ~] = size(S);
    
        phase_shifts_all = zeros(N, num_pairs);
        s2_corrected_all = zeros(N, num_pairs);
        hf = [];
    
        for p = 1:num_pairs
            s1 = S(:,p,1);
            s2 = S(:,p,2);
            [s2_corrected, phase_shifts] = ...
                CorrectBidiScanning.estimate_g_and_correct(s1, s2, num_iters, false);
            s2_corrected_all(:,p) = s2_corrected;
            phase_shifts_all(:,p) = phase_shifts;
        end
        
        % Compute consensus g
        g_consensus = median(phase_shifts_all, 2); % MEDIAN
    
        % reconstruct corrected images
        image = zeros(num_pairs*2,N);
        idxeven = mod(1:num_pairs*2,2) == 0;
        image(~idxeven,:) = S(:,:,1)';
        image(idxeven,:) = S(:,:,2)';
        image_w = image; image_w(idxeven,:) = s2_corrected_all'; % corrected line-by-line
        image_wc = CorrectBidiScanning.applycorrection(image,g_consensus); % corrected using consensus phase shift
    
        images(:,:,1) = image;
        images(:,:,2) = image_w;
        images(:,:,3) = image_wc;
        
    end

    function [s2_corrected, phase_shifts] = estimate_g_and_correct(s1, s2, num_iters, pl)
        % Parameters
        N = length(s1);
        window_sizes = round(linspace(N/4, N/20, num_iters)); % Iterative window sizes
        
        % Initialize phase shift estimate and correction result
        phase_shifts = zeros(size(s1));
        s2_corrected = s2;
        
        if pl; figure; end
        
        for iter = 1:num_iters
            win_size = window_sizes(iter);
            half_win = floor(win_size / 2);
            local_shifts = zeros(size(s1));
            
            % Compute local shifts using cross-correlation
            for i = half_win+1:N-half_win
                win_s1 = s1(i-half_win:i+half_win);
                win_s2 = s2_corrected(i-half_win:i+half_win);
                [xc, lags] = xcorr(win_s1 - mean(win_s1), win_s2 - mean(win_s2));
                [~, max_idx] = max(xc);
                local_shifts(i) = lags(max_idx);
            end
            
            % Smooth the local shifts using interpolation
            phase_shifts = phase_shifts + movmean(local_shifts,half_win);
            
            % Apply correction
            t = [1:N]';
            g_inv = t - phase_shifts; % g^-1 maps corrected time back to original
            s2_corrected = interp1(t, s2, g_inv, 'linear', 'extrap'); 
            
            if pl
                % Plot results
                subplot(num_iters, 1, iter);
                plot(t, s1, 'b', 'DisplayName', 's1'); hold on;
                plot(t, s2, 'r', 'DisplayName', 's2 (original)');
                plot(t, s2_corrected, 'g', 'DisplayName', 's2 (corrected)');
                title(['Iteration ' num2str(iter) ', Window Size: ' num2str(win_size)]);
                legend;
                hold off;
            end
        end
    end

    function img_out = applycorrection(img_in,fwd_transform)
        arguments
            img_in double
            fwd_transform double 
        end
        
        img_out = img_in;

        N = height(img_in);
        idxeven = mod(1:N,2) == 0;
        s = [img_in(idxeven,:)]'; % even lines

        % Apply correction to lines
        [N,num_lines] = size(s);
        s_corrected = zeros(N, num_lines);
        t = (1:N)';
        g_inv = t - fwd_transform; % Inverse mapping
        
        for p = 1:num_lines
            thisline = s(:, p);
            s_corrected(:, p) = interp1(t, thisline, g_inv, 'linear', 'extrap');
        end

        % reconstruct image
        img_out(idxeven,:) = s_corrected';
    end
end

methods
    function obj = CorrectBidiScanning(movie_in,reference_img,reference_meta)
            arguments
                movie_in = ''
                reference_img double = []
                reference_meta struct = struct()
            end
            obj = obj@MovieProcessing(movie_in);
            obj = obj.setReference(reference_img);
            obj.init.reference_meta = reference_meta;
            obj.init.method = 'consensus'; % this is currently the only option,
            % but it would be easy to ass a 'line-by-line' method
    end

    function obj = compute(obj, s)
        arguments
            obj 
            s logical = false
        end

        fprintf('Computing phase shifts..')

        % format data
        refimg = obj.reference_img;
        N = height(refimg);
        S(:,:,1) = [refimg(mod(1:N,2) ~= 0,:)]'; % odd lines
        S(:,:,2) = [refimg(mod(1:N,2) == 0,:)]'; % even lines
        
        % params
        num_iters = 3;
        
        % run
        [images, phase_shifts_all, g_consensus] = ...
            obj.estimate_g_and_correct_multiple(S, num_iters);

        obj.operation.images = images;
        obj.operation.phase_shifts_all = phase_shifts_all;
        obj.operation.fwd_transform_consensus = g_consensus;

        % plotting
        hf = plotOperation(obj,true);

        % optional saving
        if s; savefig('bidi_shifts.fig'); end
        
        disp(' done')
    end

    function obj = run(obj, s)
        arguments
            obj 
            s logical = false
        end

        % compute transform if necessary
        if isempty(obj.operation)
            obj = obj.compute(s);
        end

        % Signal registration start
        obj.disp_runheader()

        % initialize
        movie = obj.data_raw;
        fwd_transform = obj.operation.fwd_transform_consensus;
        movie_out = movie;
        % Preallocate output array
        [rows, cols, frames] = size(movie.stack);
        stack_out = zeros(rows, cols, frames, 'like', movie.stack);
        
        % loop frame-by-frame
        parfor i = 1:movie.nfr
            thisframe = movie.stack(:,:,i);
            stack_out(:,:,i) = obj.applycorrection(thisframe,fwd_transform);
        end
        movie_out.stack = stack_out;

        % store output
        obj.data_processed = movie_out;

        % final operations
        obj = obj.tailsequence;
    end

    function obj = setReference(obj,refimg)
        % expect refimg is double

        if isempty(refimg)
            refimg = obj.data_raw.timeavg;
        end

        obj.reference_img = refimg;
        
    end

end


end

function hf = plotOperation(obj,s)
    arguments
        obj CorrectBidiScanning
        s logical = false % optional saving
    end
    
    op = obj.operation;
    g = op.fwd_transform_consensus;
    
    hf = figure;
    subplot(231); imagesc(op.images(:,:,1)); ax1 = gca; % native img
    title('native img'); axis square
    subplot(232); imagesc(op.images(:,:,2)); ax2 = gca; % corrected line-by-line
    title('corrected line-by-line'); axis square
    subplot(233); imagesc(op.images(:,:,3)); ax3 = gca; % consensus correction
    title('consensus correction'); axis square
    linkaxes([ax1, ax2, ax3],'xy')
    subplot(234); imagesc(op.phase_shifts_all'); clim([min(g) max(g)])
    subplot(224); plot(g); axis tight
end