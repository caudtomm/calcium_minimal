function local_corr_map = computeLocalCorrelationMap(M)
    % Input: 
    %   M - 3D matrix [x, y, time] where x and y are spatial dimensions and 
    %       time is the temporal dimension.
    % Output: 
    %   local_corr_map - 2D matrix [x, y] representing the local correlation 
    %       of each pixel with its neighbors over time.

    [rows, cols, timepoints] = size(M);
    local_corr_map = zeros(rows, cols);

    % Define neighborhood (3x3 window around each pixel)
    for i = 2:rows-1
        for j = 2:cols-1
            % Extract the time series for the current pixel and its neighbors
            center_pixel_ts = squeeze(M(i, j, :));
            neighbors_ts = reshape(M(i-1:i+1, j-1:j+1, :), [], timepoints);
            
            % Compute the correlation of the center pixel with all neighbors
            correlations = corr(center_pixel_ts, neighbors_ts','rows','complete');
            
            % Take the mean correlation as the local correlation
            local_corr_map(i, j) = mean(correlations(:));
        end
    end
end
