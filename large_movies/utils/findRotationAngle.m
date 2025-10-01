function ang = findRotationAngle(v, a0)
if nargin<2, a0 = 0; end
t0 = v.CurrentTime; f = readFrame(v); try, v.CurrentTime = t0; end
if size(f,3)==3, f = rgb2gray(f); end

ang = [];
fig = figure('Name','Rotation','NumberTitle','off','MenuBar','none','ToolBar','none','Color','w');
set(fig,'Units','normalized','Position',[0.25 0.25 0.5 0.6])

ax = axes('Parent',fig,'Units','normalized','Position',[0.05 0.12 0.9 0.83]);
img = imshow(rotf(f,a0),'Parent',ax); axis(ax,'image'); ax.Visible='off';

s = uicontrol(fig,'Style','slider','Units','normalized','Position',[0.05 0.06 0.7 0.03],...
    'Min',-180,'Max',180,'Value',a0,'SliderStep',[0.1/360 1/360],'Callback',@upd);
e = uicontrol(fig,'Style','edit','Units','normalized','Position',[0.76 0.06 0.06 0.035],...
    'String',num2str(a0),'Callback',@ed);
uicontrol(fig,'Style','text','Units','normalized','Position',[0.83 0.06 0.05 0.035],...
    'String','deg','BackgroundColor','w');
b0 = uicontrol(fig,'Style','pushbutton','Units','normalized','Position',[0.05 0.02 0.08 0.035],...
    'String','0°','Callback',@(~,~) setv(0));
b9 = uicontrol(fig,'Style','pushbutton','Units','normalized','Position',[0.14 0.02 0.08 0.035],...
    'String','90°','Callback',@(~,~) setv(90));
b18= uicontrol(fig,'Style','pushbutton','Units','normalized','Position',[0.23 0.02 0.08 0.035],...
    'String','180°','Callback',@(~,~) setv(180));
b27= uicontrol(fig,'Style','pushbutton','Units','normalized','Position',[0.32 0.02 0.08 0.035],...
    'String','270°','Callback',@(~,~) setv(270));
bs = uicontrol(fig,'Style','pushbutton','Units','normalized','Position',[0.41 0.02 0.12 0.035],...
    'String','snap 90','Callback',@(~,~) setv(round(s.Value/90)*90));
cg = uicontrol(fig,'Style','checkbox','Units','normalized','Position',[0.55 0.02 0.1 0.035],...
    'String','grid','BackgroundColor','w','Callback',@gridtoggle);
ok = uicontrol(fig,'Style','pushbutton','Units','normalized','Position',[0.74 0.02 0.1 0.035],...
    'String','OK','Callback',@done);
cc = uicontrol(fig,'Style','pushbutton','Units','normalized','Position',[0.85 0.02 0.1 0.035],...
    'String','Cancel','Callback',@cancel);

set(fig,'WindowKeyPressFcn',@keys)

G = gobjects(0);
uiwait(fig)

    function upd(~,~)
        a = s.Value; e.String = num2str(a);
        img.CData = rotf(f,a);
        if ~isempty(G) && isvalid(ax), drawgrid(); end
    end
    function ed(~,~)
        a = str2double(e.String); if ~isfinite(a), a = 0; end
        a = max(-180,min(180,a)); s.Value = a; upd();
    end
    function setv(a)
        a = mod(a+180,360)-180; s.Value = a; upd();
    end
    function keys(~,ev)
        a = s.Value;
        switch ev.Key
            case 'left',  a = a - (isfield(ev,'Modifier') && any(strcmp(ev.Modifier,'shift')))*1 - (~(isfield(ev,'Modifier') && any(strcmp(ev.Modifier,'shift'))))*0.1;
            case 'right', a = a + (isfield(ev,'Modifier') && any(strcmp(ev.Modifier,'shift')))*1 + (~(isfield(ev,'Modifier') && any(strcmp(ev.Modifier,'shift'))))*0.1;
            case 'uparrow',   a = a + 1;
            case 'downarrow', a = a - 1;
            case 'r', a = 0;
            case 's', a = round(a/90)*90;
        end
        a = max(-180,min(180,a)); s.Value = a; upd();
    end
    function gridtoggle(~,~)
        if cg.Value, drawgrid(); else, delete(G(ishandle(G))); G = gobjects(0); end
    end
    function drawgrid()
        delete(G(ishandle(G))); G = gobjects(0);
        I = img.CData; [h,w] = size(I);
        d = max(32,round(min(h,w)/12));
        xs = 1:d:w; ys = 1:d:h;
        hold(ax,'on')
        for x = xs, G(end+1) = line(ax,[x x],[1 h],'Color',[0 1 0 0.25]); end %#ok<AGROW>
        for y = ys, G(end+1) = line(ax,[1 w],[y y],'Color',[0 1 0 0.25]); end %#ok<AGROW>
        hold(ax,'off')
    end
    function done(~,~)
        ang = s.Value; if isvalid(fig), uiresume(fig); delete(fig); end
    end
    function cancel(~,~)
        ang = []; if isvalid(fig), uiresume(fig); delete(fig); end
    end
end

function g = rotf(f,a)
a = mod(a,360); r = round(a/90);
if abs(a-90*r)<1e-9, g = rot90(f,r);
else, g = imrotate(f,a,'nearest','crop'); end
end
