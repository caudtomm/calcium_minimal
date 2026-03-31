classdef AnatomyAlignment < handle
    % AnatomyAlignment  Affine alignment of anatomical max-projections across subjects.
    %
    %   Loads dedicated anatomy TIFF stacks (from not_trials/*anatomy_*),
    %   generates max-intensity projections, supports interactive landmark
    %   selection, computes per-subject affine transforms, and visualises
    %   alignment quality and physiological FOV placement.
    %
    %   Typical workflow:
    %       aa = AnatomyAlignment(experiment);
    %       aa.px_sizes      = [...];   % [zoom, px_x_um, px_y_um] for anatomy
    %       aa.px_sizes_phys = [...];   % same for physiology
    %       aa.close_up_zoom = 4.0;
    %       aa.loadStacks();
    %       aa.glance();
    %       aa.selectMarkers({'OB', 'Dp', 'Vv'});
    %       aa.affineAlign();
    %       [img, px_x, px_y] = aa.showAvgAnatomy();

    properties
        experiment          Experiment          % input Experiment object
        stack               double  = []        % [H × W × N] per-subject normalised max-projections
        stack_aligned       double  = []        % [H_out × W_out × N] aligned stack, NaN background
        marker_names        cell    = {}        % {1 × P} landmark name strings
        marker_positions    double  = []        % [N × P × 2] (row, col) in original pixel space
        ref_subject         double  = []        % index into stack used as reference (empty → avg markers)
        tforms              cell    = {}        % {N × 1} affine2d transform objects (saved by affineAlign)
        px_sizes            double  = []        % [n_zoom × 3]: [zoom, px_x_um, px_y_um] for anatomy
        px_sizes_phys       double  = []        % [n_zoom × 3]: same for physiology
        close_up_zoom       double  = []        % zoom factor of the physiology recording
    end

    properties (Access = private)
        canvas struct = struct('xmin', 0, 'ymin', 0, 'width', 0, 'height', 0)
    end

    methods
        %% Constructor
        function obj = AnatomyAlignment(experiment)
            arguments
                experiment Experiment
            end
            obj.experiment = experiment;
        end

        %% loadStacks
        function obj = loadStacks(obj)
            % Load anatomy TIFFs from not_trials/*anatomy_* for each subject,
            % take max-intensity projection, normalise per subject (1st–99th pct),
            % and stack into obj.stack [H × W × N].
            N = numel(obj.experiment.traces);
            all_proj = [];

            for i = 1:N
                loc    = obj.experiment.traces{i}.subject_locations;
                nt_dir = fullfiletol(loc.subject_datapath, 'not_trials');
                files   = dir(fullfiletol(nt_dir, '*anatomy_R*'));

                if isempty(files)
                    warning('AnatomyAlignment:noFile', ...
                        'No anatomy file found for subject %d (%s)', ...
                        i, loc.subject_datapath);
                    proj = nan(height(all_proj),width(all_proj)); % assumes this is not the first file
                    all_proj(:,:,i) = proj; %#ok<AGROW>
                    continue
                end

                % Take the last file (most recent), as in scribbles.m
                fpath = fullfiletol(nt_dir, files(end).name);
                fprintf('Subject %d/%d: loading %s\n', i, N, files(end).name);

                raw  = double(Movie(fpath).stack);
                proj = max(raw, [], 3, 'omitmissing');

                % Per-subject percentile normalisation (handles outlier pixels)
                lo   = prctile(proj(:), 1);
                hi   = prctile(proj(:), 99);
                proj = (proj - lo) ./ max(hi - lo, eps);
                proj = max(0, min(1, proj));

                all_proj(:,:,i) = proj; %#ok<AGROW>
            end

            obj.stack = all_proj;
        end

        %% glance
        function glance(obj)
            % Show all max-projections packed side by side, grayscale, no labels.
            N     = size(obj.stack, 3);
            nCols = min(5, N);
            nRows = ceil(N / nCols);

            fig = figure('Color', 'k');
            set(fig, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
            tl = tiledlayout(fig, nRows, nCols, 'Padding', 'none', 'TileSpacing', 'none');

            for i = 1:N
                ax = nexttile(tl);
                imagesc(ax, obj.stack(:,:,i));
                colormap(ax, 'gray');
                clim(ax, [0 1]);
                axis(ax, 'image');
                axis(ax, 'off');
            end
        end

        %% selectMarkers
        function selectMarkers(obj, markerNames)
            % Interactively place anatomical landmarks for each subject.
            %   markerNames  – {1 × P} cell array of landmark name strings
            %
            % Controls (per subject):
            %   Left-click  : place the highlighted marker
            %   Right-click : undo the most recently placed marker
            %   Enter       : confirm all markers and move to next subject
            arguments
                obj         AnatomyAlignment
                markerNames (1,:) cell
            end

            P = numel(markerNames);
            N = size(obj.stack, 3);
            obj.marker_names     = markerNames;
            obj.marker_positions = nan(N, P, 2);

            mcolors = lines(P);

            for s = 1:N
                fig = figure('Color', 'k', ...
                    'Name', sprintf('Subject %d / %d', s, N));
                ax = axes(fig); %#ok<LAXES>
                imagesc(ax, obj.stack(:,:,s));
                colormap(ax, 'gray');
                clim(ax, [0 1]);
                axis(ax, 'image');
                axis(ax, 'off');
                hold(ax, 'on');

                h_pts  = gobjects(P, 1);
                h_text = gobjects(P, 1);

                p = 1;   % current marker index; p == P+1 → confirmation state
                while true
                    % Update window title with current state
                    if p <= P
                        set(fig, 'Name', sprintf( ...
                            'Subject %d/%d  |  Place: "%s"  [%d/%d]  |  Right-click: undo', ...
                            s, N, markerNames{p}, p, P));
                    else
                        set(fig, 'Name', sprintf( ...
                            'Subject %d/%d  |  All placed  |  Enter: confirm  |  Right-click: redo last', ...
                            s, N));
                    end

                    waitforbuttonpress;
                    key   = get(fig, 'CurrentKey');
                    stype = get(fig, 'SelectionType');

                    if p <= P && strcmp(stype, 'normal')
                        % Left-click in placement state → place current marker
                        pt = get(ax, 'CurrentPoint');
                        x  = pt(1,1);   % column
                        y  = pt(1,2);   % row
                        obj.marker_positions(s, p, :) = [y, x];

                        if isvalid(h_pts(p));  delete(h_pts(p));  end
                        if isvalid(h_text(p)); delete(h_text(p)); end
                        h_pts(p)  = plot(ax, x, y, 'o', ...
                            'Color',           mcolors(p,:), ...
                            'MarkerFaceColor', mcolors(p,:), ...
                            'MarkerSize', 9, 'LineWidth', 1.5);
                        h_text(p) = text(ax, x + 5, y - 5, markerNames{p}, ...
                            'Color',           mcolors(p,:), ...
                            'FontSize', 10, 'FontWeight', 'bold', ...
                            'Interpreter', 'none', 'BackgroundColor', 'k');
                        p = p + 1;

                    elseif strcmp(stype, 'alt')
                        % Right-click → undo most recently placed marker
                        if p == P + 1
                            p_del = P;          % undo last marker from confirm state
                        elseif p > 1
                            p_del = p - 1;      % undo previous marker during placement
                        else
                            p_del = [];         % nothing to undo at first marker
                        end

                        if ~isempty(p_del)
                            obj.marker_positions(s, p_del, :) = nan;
                            if isvalid(h_pts(p_del));  delete(h_pts(p_del));  end
                            if isvalid(h_text(p_del)); delete(h_text(p_del)); end
                            p = p_del;
                        end

                    elseif p > P && strcmp(key, 'return')
                        % Enter in confirmation state → accept and move on
                        break;
                    end
                end

                close(fig);
            end
        end

        %% affineAlign
        function obj = affineAlign(obj)
            % Compute optimal affine transforms from each subject's marker
            % positions to the reference, warp all projections onto a common
            % canvas (NaN background), and save transforms to obj.tforms and
            % aligned stack to obj.stack_aligned.

            if any(isnan(obj.marker_positions(:)))
                error('AnatomyAlignment:missingMarkers', ...
                    'Some marker positions are NaN. Complete selectMarkers first.');
            end

            [H, W, N] = size(obj.stack);

            % Reference positions [P × 2] in [col, row] = [x, y] for fitgeotrans
            if ~isempty(obj.ref_subject)
                ref_rowcol = squeeze(obj.marker_positions(obj.ref_subject, :, :));
            else
                ref_rowcol = squeeze(mean(obj.marker_positions, 1, 'omitmissing'));
            end
            ref_xy = fliplr(ref_rowcol);    % [P × 2] as [x, y]

            % Fit per-subject affine transforms
            obj.tforms = cell(N, 1);
            for i = 1:N
                mov_rowcol   = squeeze(obj.marker_positions(i, :, :));
                mov_xy       = fliplr(mov_rowcol);
                obj.tforms{i} = fitgeotrans(mov_xy, ref_xy, 'affine');
            end

            % Determine canvas size to contain every warped image entirely
            xlims = zeros(N, 2);
            ylims = zeros(N, 2);
            for i = 1:N
                [xl, yl]  = outputLimits(obj.tforms{i}, [1 W], [1 H]);
                xlims(i,:) = xl;
                ylims(i,:) = yl;
            end
            xmin  = floor(min(xlims(:,1)));
            xmax  = ceil( max(xlims(:,2)));
            ymin  = floor(min(ylims(:,1)));
            ymax  = ceil( max(ylims(:,2)));
            out_W = xmax - xmin + 1;
            out_H = ymax - ymin + 1;

            obj.canvas = struct('xmin', xmin, 'ymin', ymin, ...
                'width', out_W, 'height', out_H);

            % Common output reference (pixel size = 1, same as input)
            Rout = imref2d([out_H, out_W]);
            Rout.XWorldLimits = [xmin - 0.5,  xmin + out_W - 0.5];
            Rout.YWorldLimits = [ymin - 0.5,  ymin + out_H - 0.5];

            % Warp each projection onto the common canvas
            obj.stack_aligned = nan(out_H, out_W, N);
            for i = 1:N
                obj.stack_aligned(:,:,i) = imwarp(obj.stack(:,:,i), obj.tforms{i}, ...
                    'OutputView', Rout, 'FillValues', nan);
            end
        end

        %% showReference
        function [img, px_x_um, px_y_um] = showReference(obj)
            % Show (and return) the reference anatomy image with a 100 µm scalebar.
            %   img      – reference image matrix (no scalebar)
            %   px_x_um  – anatomy pixel size in x (µm/px)
            %   px_y_um  – anatomy pixel size in y (µm/px)
            [img, px_x_um, px_y_um] = obj.getRefImg();

            fig = figure('Color', 'k');
            ax  = axes(fig); %#ok<LAXES>
            imagesc(ax, img);
            colormap(ax, 'gray');
            axis(ax, 'image');
            axis(ax, 'off');
            c = prctile(img(:), [1 99]);
            clim(ax, c);

            addScalebar(ax, px_x_um);
        end

        %% showAvgAnatomy
        function [img, px_x_um, px_y_um] = showAvgAnatomy(obj)
            % Show (and return) the mean aligned anatomy with a 100 µm scalebar.
            %   img      – mean of stack_aligned (no scalebar, NaN-safe)
            %   px_x_um  – anatomy pixel size in x (µm/px)
            %   px_y_um  – anatomy pixel size in y (µm/px)
            img = mean(obj.stack_aligned, 3, 'omitmissing');
            [px_x_um, px_y_um] = obj.getAnatomyPxSizes();

            fig = figure('Color', 'k');
            ax  = axes(fig); %#ok<LAXES>
            imagesc(ax, img);
            colormap(ax, 'gray');
            axis(ax, 'image');
            axis(ax, 'off');
            c = prctile(img(:), [1 99]);
            clim(ax, c);

            addScalebar(ax, px_x_um);
        end

        %% showInnerAlignment
        function [img, px_x_um, px_y_um] = showInnerAlignment(obj, class_labels, cfg)
            % Overlay each subject's physiological FOV perimeter on the
            % reference anatomy, coloured by class label.
            %
            %   class_labels – (N,1) double. NaN or <1 → subject not shown.
            %                  Color = cfg.c(class_label, :).
            %   cfg          – PlotConfig (optional; default = PlotConfig())
            %
            %   img      – reference image (no scalebar, no overlays)
            %   px_x_um  – anatomy pixel size in x (µm/px)
            %   px_y_um  – anatomy pixel size in y (µm/px)
            arguments
                obj          AnatomyAlignment
                class_labels (:,1) double
                cfg          PlotConfig = PlotConfig()
            end

            [img, px_x_um, px_y_um] = obj.getRefImg();
            [px_x_anat,  px_y_anat]  = obj.getAnatomyPxSizes();

            % Physiology pixel sizes at close_up_zoom
            zoom_col = obj.px_sizes_phys(:,1);
            row = abs(zoom_col - obj.close_up_zoom) < 1e-9;
            if ~any(row)
                error('AnatomyAlignment:zoomNotFound', ...
                    'close_up_zoom %.4g not found in px_sizes_phys column 1.', ...
                    obj.close_up_zoom);
            end
            px_x_phys = obj.px_sizes_phys(row, 2);
            px_y_phys = obj.px_sizes_phys(row, 3);

            % Physiology FOV size in anatomy-image pixels
            % (512 x 256 physiology px, physically square in µm)
            fov_w_anatpx = 512 * px_x_phys / px_x_anat;
            fov_h_anatpx = 256 * px_y_phys / px_y_anat;

            [H_raw, W_raw, ~] = size(obj.stack);

            % Display reference on canvas
            fig = figure('Color', 'k');
            ax  = axes(fig); %#ok<LAXES>
            imagesc(ax, img);
            colormap(ax, 'gray');
            axis(ax, 'image');
            axis(ax, 'off');
            c = prctile(img(:), [1 99]);
            clim(ax, c);
            hold(ax, 'on');

            N = numel(class_labels);
            for i = 1:N
                lbl = class_labels(i);
                if isnan(lbl) || lbl < 1; continue; end

                % FOV corners centred on subject i's raw anatomy (world coords, [x y])
                cx = W_raw / 2;
                cy = H_raw / 2;
                hw = fov_w_anatpx / 2;
                hh = fov_h_anatpx / 2;
                corners_xy = [cx-hw, cy-hh;     % top-left
                              cx+hw, cy-hh;     % top-right
                              cx+hw, cy+hh;     % bottom-right
                              cx-hw, cy+hh;     % bottom-left
                              cx-hw, cy-hh];    % close perimeter

                % Transform corners into canvas world coordinates
                tc = transformPointsForward(obj.tforms{i}, corners_xy);

                % Convert world → canvas pixel indices
                tc_px_x = tc(:,1) - obj.canvas.xmin + 1;
                tc_px_y = tc(:,2) - obj.canvas.ymin + 1;

                color = cfg.c(round(lbl), :);
                plot(ax, tc_px_x, tc_px_y, '-', 'Color', color, 'LineWidth', 1.5);
            end

            addScalebar(ax, px_x_um);
        end
    end

    %% Private helpers
    methods (Access = private)
        function [img, px_x_um, px_y_um] = getRefImg(obj)
            [px_x_um, px_y_um] = obj.getAnatomyPxSizes();
            if ~isempty(obj.ref_subject)
                img = obj.stack_aligned(:,:,obj.ref_subject);
            else
                img = mean(obj.stack_aligned, 3, 'omitmissing');
            end
        end

        function [px_x_um, px_y_um] = getAnatomyPxSizes(obj)
            row = abs(obj.px_sizes(:,1) - 1.0) < 1e-9;
            if ~any(row)
                error('AnatomyAlignment:noZoom1', ...
                    'Zoom 1.0 entry not found in px_sizes column 1.');
            end
            px_x_um = obj.px_sizes(row, 2);
            px_y_um = obj.px_sizes(row, 3);
        end
    end
end

%% ── File-scope helpers ───────────────────────────────────────────────────────

function addScalebar(ax, px_x_um)
    % Overlay a 100 µm white scalebar at the bottom-right of ax.
    scalebar_px = 100 / px_x_um;

    xl      = xlim(ax);
    yl      = ylim(ax);
    x_range = xl(2) - xl(1);
    y_range = yl(2) - yl(1);

    x2    = xl(1) + 0.95 * x_range;
    x1    = x2 - scalebar_px;
    y_bar = yl(1) + 0.92 * y_range;    % near bottom (y increases downward in imagesc)
    y_txt = yl(1) + 0.87 * y_range;

    line(ax, [x1, x2], [y_bar, y_bar], 'Color', 'w', 'LineWidth', 3);
    text(ax, (x1+x2)/2, y_txt, '100 µm', ...
        'Color', 'w', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 10, 'FontWeight', 'bold');
end
