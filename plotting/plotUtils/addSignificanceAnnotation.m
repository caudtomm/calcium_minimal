function h = addSignificanceAnnotation(ax, x, y, p, opts)
    % addSignificanceAnnotation - Add significance stars or brackets to a plot
    %
    % Usage:
    %   addSignificanceAnnotation(ax, x, y, p)
    %   addSignificanceAnnotation(ax, [x1 x2], y, p, 'style', 'bracket')
    %
    % Inputs:
    %   ax - axes handle (or [] to use gca)
    %   x  - x position (scalar) or [x1 x2] for bracket
    %   y  - y position for annotation
    %   p  - p-value (will be converted to stars)
    %
    % Options:
    %   style    - 'stars' (default) or 'bracket'
    %   color    - text/line color (default: 'k')
    %   fontSize - font size for stars (default: 12)
    %   offset   - vertical offset for bracket arms (default: 0.02 * y-range)
    %
    % Output:
    %   h - handle(s) to created graphics objects
    %
    % Example:
    %   % Simple stars above a point
    %   addSignificanceAnnotation(gca, 2, 0.5, 0.003);
    %
    %   % Bracket between two groups
    %   addSignificanceAnnotation(gca, [1 2], 0.6, 0.01, 'style', 'bracket');

    arguments
        ax
        x double
        y double
        p double
        opts.style char {mustBeMember(opts.style, {'stars', 'bracket'})} = 'stars'
        opts.color = 'k'
        opts.fontSize double = 12
        opts.offset double = []
    end

    if isempty(ax)
        ax = gca;
    end

    % Convert p-value to stars
    stars = statsUtils.pvalToStars(p);

    % Don't annotate if not significant
    if strcmp(stars, 'n.s.') || isempty(stars)
        h = [];
        return;
    end

    % Calculate offset if not provided
    if isempty(opts.offset)
        yl = ylim(ax);
        opts.offset = 0.02 * (yl(2) - yl(1));
    end

    hold(ax, 'on');

    if strcmp(opts.style, 'stars') || isscalar(x)
        % Simple stars annotation
        x_pos = mean(x);  % center if x is a range
        h = text(ax, x_pos, y, stars, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'FontSize', opts.fontSize, ...
            'Color', opts.color);
    else
        % Bracket style: line with stars above
        x1 = x(1);
        x2 = x(2);
        x_mid = (x1 + x2) / 2;
        y_bracket = y;
        y_arms = y - opts.offset;

        % Draw bracket
        h_line = plot(ax, [x1 x1 x2 x2], [y_arms y_bracket y_bracket y_arms], ...
            'Color', opts.color, 'LineWidth', 1);

        % Add stars
        h_text = text(ax, x_mid, y_bracket + opts.offset/2, stars, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'FontSize', opts.fontSize, ...
            'Color', opts.color);

        h = [h_line; h_text];
    end
end
