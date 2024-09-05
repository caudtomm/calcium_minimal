function ROImap = imageSequenceGUI(images,localCorr,ROImap,GUItitle)

    if nargin<2; error('Not enough input arguments.'); end

    % Reformat images if necessary
    if isnumeric(images); images = format2cell(images); end
    if isnumeric(localCorr); localCorr = format2cell(localCorr); end

    % set default GUI title
    if ~exist("GUItitle",'var') || isempty(GUItitle)
        GUItitle = 'Image Sequence Viewer';
    end

    numImages = numel(images);

    if ~exist("ROImap",'var') || isempty(ROImap)
        % Initialize ROImap with zeros, same size as the first image
        ROImap = zeros(size(images{1}, 1), size(images{1}, 2));
    end
    

    % Normalize images for visualization
    images = normalizeImages(images);
    localCorr = normalizeImages(localCorr);

    %% Generate GUI

    % Create the main figure
    fig_h = 600;
    fig_w = 800;
    fig = figure('Name', GUItitle, 'NumberTitle', 'off', ...
                 'MenuBar', 'none', 'ToolBar', 'none', 'Position', [100, 100, fig_w, fig_h], ...
                 'KeyPressFcn', @keyPressCallback, 'WindowScrollWheelFcn', @mouseScrollCallback, ...
                 'WindowButtonDownFcn', @detectMouseClick);
    fig_szfactor = fig_w/fig_h;

    % Main image axes
    axMain = axes('Parent', fig, 'Units', 'normalized', ...
                  'Position', [0.2, 0.05, 0.75, 0.9]);

    % Initialize the main image
    currentImageIndex = 1;
    hMainImage = imshow(images{currentImageIndex}, 'DisplayRange', [0 1], 'Parent', axMain);
    axis(axMain,'square')
    hold(axMain, 'on'); % Allow overlaying drawings
    fullXLim = get(axMain, 'XLim'); originalXLim = fullXLim;
    fullYLim = get(axMain, 'YLim'); originalYLim = fullYLim;
    visualMode = 'raw_roi';


    % Panel for thumbnails
    thumbPanel_h = 0.9;
    thumbPanel_w = 0.15;
    thumbPanel = uipanel('Parent', fig, 'Units', 'normalized', ...
                         'Position', [0.05, 0.05, thumbPanel_w, thumbPanel_h]);
    thumbPfactor = thumbPanel_w/thumbPanel_h;
    szfactor = thumbPfactor/fig_szfactor;

    % Create thumbnails in the panel
    thumbAxes = gobjects(1, numImages);
    thumbImages = gobjects(1, numImages);
    thumbRects = gobjects(1, numImages); % Array to hold rectangles for highlighting
    thumbSize = .85; % Size of each thumbnail
    initialOffset = .05; % Initial offset for scrolling
    offset = initialOffset;
    for i = 1:numImages
        thumbAxes(i) = axes('Parent', thumbPanel, 'Units', 'normalized', ...
                            'Position', [0.05, thumbPfactor*(i-1) + offset, thumbSize, thumbSize*thumbPfactor], ...
                            'ButtonDownFcn', @(~,~) updateMainImage(i));
        thumbImages(i) = imshow(images{i}, 'DisplayRange', [0 1], 'Parent', thumbAxes(i));
        set(thumbImages(i), 'ButtonDownFcn', @(~,~) updateMainImage(i));
        axis square

        % Overlay the index number on the thumbnail
        [thumbHeight,thumbWidth] = size(thumbImages(i).CData);
        ht = text(thumbAxes(i), thumbWidth/2, thumbHeight/2, num2str(i), 'Color', 'yellow', ...
             'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

        % Create a rectangle around each thumbnail (initially invisible)
        thumbRects(i) = rectangle(thumbAxes(i), 'Position', [0, 0, 1, 1], ...
                                  'EdgeColor', 'none', 'FaceColor', 'none', 'LineWidth', 20);
    end
    highlightThumbnail(currentImageIndex);

    updateMainImage(currentImageIndex);

    freehandROI = [];    % Handle for the freehand drawing

    % Pause the execution until the GUI is closed
    waitfor(fig);

    %% functions and callbacks

    % reformat double array input to cell array (numel = n images)
    function outcell = format2cell(indouble)
        indouble(indouble==-1)=nan;
        outcell = {};
        for i_im = 1:size(indouble,3)
            outcell{end+1} = indouble(:,:,i_im);
        end
    end

    % normalize images for visualization
    function outcell = normalizeImages(incell)
        outcell = cell(size(incell));
        for i_im = 1:numel(incell)
            thisim = incell{i_im};
            
            % equalize histogram
            thisim = adapthisteq(mat2gray(thisim));
            
            % normalize
            normfactor = quantile(thisim(:),.99);
            minval = nanmin(thisim,[],'all');
            thisim = (thisim - minval) ./ normfactor;
            
            outcell{i_im} = thisim;
        end
    end

    % Callback to update the main image when a thumbnail is clicked
    function updateMainImage(index)
        currentImageIndex = index;

        switch visualMode
            case 'raw'
                hMainImage = imshow(images{currentImageIndex}, ...
                    'DisplayRange', [0 1], 'Parent', axMain, 'Colormap', gray);
                hold(axMain, 'on'); % Ensure drawings are preserved
            case 'raw_roi'
                hMainImage = imshow(images{currentImageIndex}, ...
                    'DisplayRange', [0 1], 'Parent', axMain, 'Colormap', gray);
                hold(axMain, 'on'); % Ensure drawings are preserved
                overlayROIs();
            case 'localcorr'
                hMainImage = imshow(localCorr{currentImageIndex}, ...
                    'DisplayRange', [0 1], 'Parent', axMain, 'Colormap', gray);
                hold(axMain, 'on'); % Ensure drawings are preserved
            case 'localcorr_roi'
                hMainImage = imshow(localCorr{currentImageIndex}, ...
                    'DisplayRange', [0 1], 'Parent', axMain, 'Colormap', gray);
                hold(axMain, 'on'); % Ensure drawings are preserved
                overlayROIs();
            case 'roimap'
                hMainImage = imshow(ROImap/max(ROImap,[],'all'), ...
                    'Parent', axMain, 'Colormap', parula);
                hold(axMain, 'on'); % Ensure drawings are preserved
        end
        axis(axMain,'square')
        highlightThumbnail(currentImageIndex);
        

        % Store the new original axis limits for the new image
        originalXLim = get(axMain, 'XLim');
        originalYLim = get(axMain, 'YLim');
    end

    % Function to handle mouse scroll
    function mouseScrollCallback(~, event)
        mousePosition = get(fig, 'CurrentPoint');  % Get the current mouse position

        % Determine if the mouse is over the thumbnail panel
        if mousePosition(1) >= 0.05 * fig.Position(3) && mousePosition(1) <= 0.25 * fig.Position(3) && ...
           mousePosition(2) >= 0.1 * fig.Position(4) && mousePosition(2) <= 0.9 * fig.Position(4)
            % If over the thumbnail panel, scroll the thumbnails
            mouseScrollThumbnails(event);
        else
            % If over the main image, zoom in/out
            mouseScrollZoom(event);
        end
    end

    % Function to handle mouse scroll for thumbnails
    function mouseScrollThumbnails(event)
        % Adjust the offset based on scroll direction
        scrollDirection = event.VerticalScrollCount;
        scrollSpeed = szfactor;
        shiftThumbs(offset + scrollDirection * scrollSpeed);
    end

    % Function to handle mouse scroll for zooming in/out on the main image
    function mouseScrollZoom(event)
        % Get the current axis limits
        axLimits = axis(axMain);
        % Get the cursor position within the axis
        cursorPoint = get(axMain, 'CurrentPoint');
        xCursor = cursorPoint(1, 1);
        yCursor = cursorPoint(1, 2);

        % Determine the zoom factor
        zoomFactor = 1.2;
        if event.VerticalScrollCount > 0
            % Zoom out
            newXLim = [xCursor - (xCursor - axLimits(1)) * zoomFactor, ...
                       xCursor + (axLimits(2) - xCursor) * zoomFactor];
            newYLim = [yCursor - (yCursor - axLimits(3)) * zoomFactor, ...
                       yCursor + (axLimits(4) - yCursor) * zoomFactor];
        else
            % Zoom in
            newXLim = [xCursor - (xCursor - axLimits(1)) / zoomFactor, ...
                       xCursor + (axLimits(2) - xCursor) / zoomFactor];
            newYLim = [yCursor - (yCursor - axLimits(3)) / zoomFactor, ...
                       yCursor + (axLimits(4) - yCursor) / zoomFactor];
        end

        % Set the new axis limits
        axis(axMain, [newXLim, newYLim]);
    end
   
   % Function to detect mouse clicks
    function detectMouseClick(~, ~)
        clickType = get(fig, 'SelectionType');
        mousePosition = get(fig, 'CurrentPoint');
        axPosition = getpixelposition(axMain, true);

        % Check if the click is within the main image axes
        if mousePosition(1) >= axPosition(1) && mousePosition(1) <= axPosition(1) + axPosition(3) && ...
           mousePosition(2) >= axPosition(2) && mousePosition(2) <= axPosition(2) + axPosition(4)
            if strcmp(clickType, 'alt')  % Right-click to reset zoom
                axis(axMain, [fullXLim, fullYLim]);
            end
        end
    end

    % Function to start freehand drawing
    function startFreehandDrawing(~, ~)
        % Create freehand ROI
        freehandROI = drawfreehand(axMain, 'Color', 'r', 'LineWidth', 2);

        if ~isempty(freehandROI) && isvalid(freehandROI)
            % Get the mask of the drawn ROI
            mask = createMask(freehandROI, hMainImage);

            % Find the next positive integer not used in ROImap
            nextLabel = max(ROImap(:)) + 1;

            % Update ROImap with the new ROI
            ROImap(mask) = nextLabel;

            % Optionally, overlay the ROI on the image for visualization
            visboundaries(axMain, mask, 'Color', 'g', 'LineWidth', 1);

            % Delete the freehand ROI object
            delete(freehandROI);
            freehandROI = [];
        end
    end

    function deleteROI()    
        % Notify the user that they can click on an ROI to delete it
        disp('Click on an ROI to delete it. Press ESC to exit deletion mode without deleting.');
    
        % Wait for the user to click on the image to select an ROI
        roiToDelete = [];
        [x, y] = ginput(1); % Get a single mouse click position
    
        % Convert the clicked point to a binary mask
        if ~isempty(x) && ~isempty(y)
            % Get the ROI label at the clicked position
            clickedLabel = ROImap(round(y), round(x));
    
            if clickedLabel > 0
                % Create a mask of the ROI to delete
                roiToDelete = (ROImap == clickedLabel);
    
                % Remove the ROI from the ROImap
                ROImap(roiToDelete) = 0;
    
                % Update the display by redrawing the image without the
                % deleted ROI
                updateMainImage(currentImageIndex);
            else
                disp('No ROI found at the clicked position.');
            end
        end
    end

    function overlayROIs()
        % Overlay remaining ROIs
        remainingLabels = unique(ROImap(ROImap > 0));
        for i_r = 1:length(remainingLabels)
            mask = (ROImap == remainingLabels(i_r));
            visboundaries(axMain, mask, 'Color', 'g', 'LineWidth', 1);
        end
    end

    % Function to highlight the selected thumbnail
    function highlightThumbnail(index)
        for j = 1:numImages
            if j == index
                set(thumbRects(j), 'EdgeColor', 'red'); % Highlight the selected thumbnail
            else
                set(thumbRects(j), 'EdgeColor', 'none'); % Remove highlight from others
            end
        end
    end

    % Update the position of each thumbnail based on the new offset
    function shiftThumbs(newoffset)
        maxUp = -szfactor/thumbSize* numImages + .1;
        maxDown = initialOffset;

        offset = newoffset;
        % Limit the scrolling to the image sequence boundaries
        if offset > maxDown
            offset = maxDown;
        elseif offset < maxUp
            offset = maxUp;
        end

        for j = 1:numImages
            set(thumbAxes(j), 'Position', [0.05, thumbPfactor*(j-1) + offset, thumbSize, thumbSize*thumbPfactor]);
        end
    end

    % Function to handle key press events (up/down arrow keys)
    function keyPressCallback(~, event)
        switch event.Key
            case 'downarrow'
                if currentImageIndex > 1
                    currentImageIndex = currentImageIndex - 1;
                    updateMainImage(currentImageIndex);
                    shiftThumbs(offset + thumbPfactor);
                end
            case 'uparrow'
                if currentImageIndex < numImages
                    currentImageIndex = currentImageIndex + 1;
                    updateMainImage(currentImageIndex);
                    shiftThumbs(offset - thumbPfactor);
                end
            case 'f'  % Enter drawing mode when 'F' is pressed
                disp('Entered drawing mode: Draw your ROI by left-clicking and dragging. Release to finalize.');
                startFreehandDrawing();
            case 'd'  % Enter deletion mode when 'D' is pressed
                disp('Entered deletion mode: Delete an ROI by clicking on it.');
                deleteROI();
            case 'q'
                visualMode = 'raw';
                updateMainImage(currentImageIndex);
            case 'w'
                visualMode = 'raw_roi';
                updateMainImage(currentImageIndex);
            case 'e'
                visualMode = 'localcorr';
                updateMainImage(currentImageIndex);
            case 'r'
                visualMode = 'localcorr_roi';
                updateMainImage(currentImageIndex);
            case 't'
                visualMode = 'roimap';
                updateMainImage(currentImageIndex);
        end
    end

end
