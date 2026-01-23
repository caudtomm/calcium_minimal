function y = plotLineNShade(t, mu, err, c, cfg)
    t = t(:)';
    mu = mu(:);
    err = err(:);
    
    fill([t fliplr(t)], [mu - err; flipud(mu + err)]', ...
        c, 'FaceAlpha', 0.4, 'EdgeColor', 'none');
    hold on
    y = plot(t,mu,'Color',c,'LineStyle',cfg.lineStyle,'LineWidth',cfg.lineWidth);
end
