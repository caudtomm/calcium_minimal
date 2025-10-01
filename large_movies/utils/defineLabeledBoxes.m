function m = defineLabeledBoxes(v, ang)
% defineLabeledBoxes - draw & name ROIs on a rotated preview
% INPUT
%   v   : VideoReader
%   ang : rotation angle (deg, CCW). Default 0
% OUTPUT
%   m   : containers.Map name -> struct('Xs',[xmin xmax], 'Ys',[ymin ymax])

if nargin<2, ang = 0; end

% preview frame (robust to EOF)
t0 = v.CurrentTime;
f = try_read_frame(v);
v.CurrentTime = t0;
if size(f,3)==3, f = rgb2gray(f); end
f = rot_frame(f, ang);

% figure + axes
hf = figure('Name','defineLabeledBoxes','NumberTitle','off');
ax = axes('Parent',hf); imshow(f,'Parent',ax); axis(ax,'image'); ax.Visible='off';
title(ax, sprintf('Draw ROI (ang=%.3f°). Double-click, then name', ang))

[H,W] = size(f);
m = containers.Map('KeyType','char','ValueType','any');
idx = 0;

while ishghandle(hf)
    % draw ROI
    h = imrect(ax); p = wait(h);
    if isempty(p) || ~ishghandle(hf), break; end

    % clamp & convert
    x1 = max(1, min(W, p(1)));  y1 = max(1, min(H, p(2)));
    x2 = max(1, min(W, p(1)+p(3)));  y2 = max(1, min(H, p(2)+p(4)));
    if ~(x2>x1 && y2>y1), uiwait(errordlg('ROI must have positive size','ROI')); continue; end

    % prompt for name (reprompt until valid)
    def = sprintf('roi_%02d', idx+1);
    nm = prompt_unique_name(m, def);
    if isempty(nm), break; end
    idx = idx + 1;

    % store
    m(nm) = struct('Xs',[x1 x2], 'Ys',[y1 y2]);

    % overlay
    rectangle('Position',[x1 y1 (x2-x1) (y2-y1)],'EdgeColor','g','LineWidth',1.2,'Parent',ax)
    text(x1, max(1,y1-5), nm, 'Color','g','FontWeight','bold','Interpreter','none','Parent',ax)

    % continue?
    q = questdlg('Add another ROI?','ROIs','Yes','No','Yes');
    if ~strcmp(q,'Yes'), break; end
end

if ishghandle(hf), close(hf); end
end

function f = try_read_frame(v)
% best-effort single frame
try
    f = readFrame(v);
catch
    v.CurrentTime = 0;
    f = readFrame(v);
end
end

function nm = prompt_unique_name(m, def)
nm = '';
while true
    a = inputdlg('ROI name (unique, folder-safe):','Name',[1 50],{def});
    if isempty(a), nm = ''; return; end
    nm = regexprep(strtrim(a{1}),'[^A-Za-z0-9_\-]','_');
    if isempty(nm)
        uiwait(errordlg('Name cannot be empty','Name')); continue
    end
    if isKey(m,nm)
        uiwait(errordlg('Name already exists','Name')); continue
    end
    if any(strcmp(nm, {'.','..'}))
        uiwait(errordlg('Invalid folder name','Name')); continue
    end
    break
end
end

function g = rot_frame(f, ang)
a = mod(ang,360);
r = round(a/90);
if abs(a-90*r)<1e-9
    g = rot90(f, r);
else
    g = imrotate(f, ang, 'nearest', 'crop'); % fast + constant size
end
end
