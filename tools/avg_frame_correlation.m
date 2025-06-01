function corr_curve = avg_frame_correlation(data)
    [rows, cols, frames] = size(data);
    corr_curve = zeros(frames - 1, 1);
    
    for i = 1:frames-1
        frame1 = reshape(data(:,:,i), [], 1);
        frame2 = reshape(data(:,:,i+1), [], 1);
        corr_curve(i) = corr(frame1, frame2);
    end
end