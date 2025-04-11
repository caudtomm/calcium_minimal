function [sorted_images, idx, results, hf] = sortImgsPOW(images, pl)
% output: low->high freq power bias
arguments
    images cell % containing double images
    pl logical = false % optional plotting
end

% knobs
th = .1; % normalized power thresh -> for bias score

%init vars
bias_scores = zeros(numel(images),1);
[Ps, rpps] = deal(cell(numel(images),1));

for i = 1:numel(images)
    img = images{i};

    % Handle NaNs
    if any(isnan(img(:)))
        m = mean(img,'all','omitmissing');
        img(isnan(img)) = m;
    end
    
    % FFT
    F = fft2(img);
    F = fftshift(F);
    P = abs(F).^2;

    % Radial frequency map
    [h, w] = size(P);
    [X, Y] = meshgrid(1:w, 1:h);
    cx = (w+1)/2; cy = (h+1)/2;
    R = sqrt((X - cx).^2 + (Y - cy).^2);

    % Radial power profile
    Rint = round(R);
    rpp = accumarray(Rint(:)+1, P(:)) ./ accumarray(Rint(:)+1, 1);
    rpp = rpp(2:end); % first value is usually NaN or otherwise unreliable

    % Normalize radius to [0,1]
    R = R ./ max(R(:));

    % Define low/high threshold (e.g. 0.5)
    low = P(R <= th);
    high = P(R > th);
    
    % Compute bias (low power ratio)
    bias_scores(i) = sum(low(:)) / (sum(low(:)) + sum(high(:)));

    % export spectra (for plotting)
    Ps{i} = P;
    rpps{i} = rpp;
end

% Sort images by low-frequency power bias (descending)
[~, idx] = sort(bias_scores, 'descend');
sorted_images = images(idx);

% extract additional metrics
target_len = max(cellfun(@numel, rpps));
[M,rpps_samesize,rppnorms] = deal(zeros(numel(images),target_len));
for i = 1:numel(images)
    rpp = rpps{idx(i)};
    n = length(rpp);
    x = linspace(1,target_len,n);
    M(i,:) = interp1(x,log(rpp),1:target_len,'linear','extrap');
    rpps_samesize(i,:) = interp1(x,rpp,1:target_len,'linear','extrap');
end
[centroids, slopes, KLs] = deal(zeros(numel(images),1));
rpps_samesize(rpps_samesize<0) = 0;
ref = mean(rpps_samesize,1,'omitmissing'); ref = ref / sum(ref,"omitmissing");
for i = 1:numel(images)
    rpp = rpps{idx(i)};
    n = length(rpp);
    
    % calc centroids
    r = (0:n-1)';
    centroids(i) = sum(r.*log(rpp),'omitmissing')/sum(log(rpp),'omitmissing');

    % calc spectral slope (in loglog space)
    r = (1:n)';
    valid = rpp>0;
    p = polyfit(log(r(valid)), log(rpp(valid)),1);
    slopes(i) = p(1);

    % calc KL divergence from reference spectrum
    % normalize
    rpp = rpps_samesize(i,:);
    rppnorm = rpp / sum(rpp,'omitmissing');
    KLs(i) = sum(rppnorm .* log(rppnorm ./ ref), 'omitmissing');
    rppnorms(i,:) = rppnorm;
end

% return results
results.bias_scores = bias_scores;
results.pow_spectra = Ps;
results.mesh_normalized = R;
results.rad_pow_profiles = rpps;
results.interp_spectra.raw = rpps_samesize;
results.interp_spectra.log = M;
results.interp_spectra.norm = rppnorms;
results.centroids = centroids;
results.loglogslopes = slopes;
results.normref = ref;
results.KLdivergence = KLs;

% PLOTTING
hf = gobjects(1,4);
if ~pl; return; end

hf(1) = figure; % power spectra
ncols = floor(numel(images)/2)+1;
nrows = floor(numel(images)/(ncols-1));

t = tiledlayout(nrows, ncols, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:numel(images)
    nexttile
    imagesc(log(1+Ps{idx(i)})); axis image off
    colormap('jet'); colorbar
    title(['Log PowSp (#',num2str(idx(i)),')'])
end

hf(2) = figure; % distribution of scores
x = 1+ .5*(rand(numel(bias_scores),1)-.5);
scatter(x,bias_scores,"filled");
ylim([0 1]); xlim([.5 1.5])
xticks(1); xticklabels({'power bias scores'});
for i = 1:numel(bias_scores) % Add index labels
    text(x(i) + 0.03, bias_scores(i), num2str(idx(i)), 'FontSize', 8,'Color','b');
end
set(hf(2),"Position",[1 1 200 300])

hf(3) = figure; % radial power profiles (normalized) and spectral centroid
subplot(141)
n = target_len;
x = 0:1/n:1-1/n;
semilogy(x,rpps_samesize)
xlabel('Radius')
ylabel('Average Power')
title('Radial Power Profiles')

subplot(142)
n = target_len;
x = 0:1/n:1-1/n;
semilogy(x,rppnorms); hold on
b = semilogy(x,ref,'r--','LineWidth',1.5);
legend(b,'reference');
xlabel('Radius')
ylabel('Average Power')
title('Norm. Profiles')


subplot(143) % plot again with centroids
imagesc(M); colorbar % spectral profiles
hold on
y = 1:numel(images);
b = scatter(centroids, y, 50, 'r' , 'filled'); % centroids
legend(b,'centroids')
xticks([])
xlabel('Radius')
yticklabels(idx)
ylabel('Image #')
title('Log Radial Power Profile')

subplot(144) % spectral slope distribution
x = 1+ .5*(rand(numel(slopes),1)-.5);
scatter(x,slopes,"filled");
xlim([.5 1.5])
xticks([]); title('Spectral Slopes')
for i = 1:numel(slopes) % Add index labels
    text(x(i) + 0.03, slopes(i), num2str(idx(i)), 'FontSize', 8,'Color','b');
end
set(hf(3),"Position",[1 1 800 500])

hf(4) = figure; % comparison between metrics
subplot(131)% bias scores vs spectral slope
scatter(1:numel(images),slopes, 'filled'); axis square
xticks(1:numel(images)); xticklabels(idx);
xlabel('Decreasing bias score')
ylabel('Spectral Slope')
hold on
rsquared = fitlm(1:numel(images),slopes).Rsquared.Ordinary;
text(1,mean(slopes),['R^2 = ',num2str(rsquared)],'FontSize',8,'Color','r');

subplot(132) % bias scores vs centroids
scatter(1:numel(images),centroids, 'filled'); axis square
xticks(1:numel(images)); xticklabels(idx);
xlabel('Decreasing bias score')
ylabel('Spectral Centroids')
hold on
rsquared = fitlm(1:numel(images),centroids).Rsquared.Ordinary;
text(1,mean(centroids),['R^2 = ',num2str(rsquared)],'FontSize',8,'Color','r');

subplot(133) % bias scores vs KL divergence
scatter(1:numel(images),KLs, 'filled'); axis square
xticks(1:numel(images)); xticklabels(idx);
xlabel('Decreasing bias score')
ylabel('KL Divergence')
hold on
rsquared = fitlm(1:numel(images),KLs).Rsquared.Ordinary;
text(1,mean(KLs),['R^2 = ',num2str(rsquared)],'FontSize',8,'Color','r');

set(hf(4),'Position',[1 1 800 300])
end
