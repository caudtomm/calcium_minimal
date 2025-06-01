function save_figures_to_pdf(hf, filename)
%SAVE_FIGURES_TO_PDF Save multiple figures as pages in a single PDF
%   hf: array of figure handles
%   filename: name of the output PDF file (with or without '.pdf' extension)

if nargin < 2 || isempty(filename)
    filename = 'figures.pdf';
end
[~,~,ext] = fileparts(filename);
if isempty(ext)
    filename = [filename '.pdf'];
end

% Loop through figure handles and append each one to the PDF
for i = 1:numel(hf)
    if isempty(hf(i)) || ~isgraphics(hf(i)); continue; end
    exportgraphics(hf(i), filename, 'Append', i > 1);
end

disp(['Saved all figures to ' filename]);
end