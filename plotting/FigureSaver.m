classdef FigureSaver < handle
    properties
        outputfolder char = 'FigureOutput'
        outputfile char = 'output.pdf'
        title char = ''
        config PlotConfig = PlotConfig() % Optional PlotConfig object
    end

    properties (Access = private)
        pageCounter double = 0
    end

    methods
        function append(obj, hf)
            if ~exist(obj.outputfolder, 'dir')
                mkdir(obj.outputfolder);
            end

            if ~all(ishandle(hf)) || ~all(strcmp(get(hf, 'Type'), 'figure'))
                error('Input must be a figure handle or array of figure handles.');
            end

            for i = 1:numel(hf)
                if ~isvalid(hf(i)), continue; end
                figure(hf(i));  % Bring to foreground
                obj.pageCounter = obj.pageCounter + 1;

                % Construct figure identifier and title
                figNumStr = sprintf('Figure_%03d', obj.pageCounter);
                headerText = sprintf('%s: %s', figNumStr, obj.title);

                % Construct full file paths
                basePath = fullfile(obj.outputfolder, figNumStr);
                svgPath = [basePath, '.svg'];
                figPath = [basePath, '.fig'];
                pdfPath = fullfile(obj.outputfolder, obj.outputfile);

                % Save as individual files (without title overlay)
                saveas(hf(i), figPath);
                saveas(hf(i), svgPath);

                % Increase figure height slightly to make space for the
                % annotation
                ax = findall(hf(i), 'Type', 'axes');
                for k = 1:numel(ax)
                    set(ax(k), 'Units', 'pixels'); % Prevent scaling
                end
                pos = get(hf(i), 'Position');
                pos(4) = pos(4) + 50;  % increase height
                hf(i).Position = pos;

                % Add annotation title (for PDF only)
                txtColor = obj.config.textcol;
                annotation(hf(i), 'textbox', [0 0.95 1 0.05], ...
                    'String', headerText, ...
                    'HorizontalAlignment', 'center', ...
                    'FontWeight', 'bold', ...
                    'EdgeColor', 'none', ...
                    'Color', txtColor, ...
                    'FontSize', 12);

                % Append to PDF
                try
                    exportgraphics(hf(i), pdfPath, ...
                        'ContentType', 'vector', ...
                        'BackgroundColor',obj.config.bgcol, ...
                        'Append', true);
                catch ME
                    warning('Failed to export Figure %d to PDF: %s', obj.pageCounter, ME.message);
                end
            end
        end
    end
end
