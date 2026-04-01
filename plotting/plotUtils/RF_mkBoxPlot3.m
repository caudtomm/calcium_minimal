function RF_mkBoxPlot3(datacells,colormatrix,xvals,boxwidth,sdlength,linw,meanmarksize,boxfillcolormat);

%datacells: cell array of size 1 x n. Each cell contains one column vector,
%which contains a set of individual datapoints that are then averaged/medianed.
%Successive cells are successive points on the xaxis. E.g., each cell may
%contain a vector with activity value at a time point tx, and successive
%cells represent successive time points. Boxes will then show mean/median
%activity across neurons as a function of time.
%
%To get started you can try:
%RF_mkBoxPlot3(datacells,[],[],[],1,1,5,[]);
%

if nargin<2;
    colormatrix=[];
end

if nargin<3;
    xvals=[];
    boxwidth=[];
    sdlength=[];
    linw=1;
    meanmarksize=0;
    boxfillcolormat=[];
end
if nargin<4;
    boxwidth=[];
    sdlength=[];
    linw=1;
    meanmarksize=0;
    boxfillcolormat=[];
end
if nargin<5;
    sdlength=[];
    linw=1;
    meanmarksize=0;
    boxfillcolormat=[];
end

if isempty(meanmarksize);
    meanmarksize=0;
end

boxwdefault=0.8;
percentiledefs=[25 75];
linw1=linw;

if ~iscell(datacells);
    help={};
    for a1=1:size(datacells,2);
        help{a1}=squeeze(datacells(:,a1));
    end
    datacells=help;
end

if isempty(xvals);
    xvals=1:length(datacells);
%     xvals=1:size(datacells,2);
end

if isempty(colormatrix);
    colormatrix=[0 0 0];
end
if isempty(boxfillcolormat);
    boxfillcolormat=[1 1 1];
    nofacecolor=1;
else
    nofacecolor=0;
end

if isempty(sdlength);
    sdlength=0;
end


ngroups=length(datacells);
for a1=1:ngroups;
    m=mean(datacells{a1});
    sd=std(datacells{a1});
    med=median(datacells{a1});
    pctles=prctile(datacells{a1},percentiledefs);
    if ngroups==1;
        if isempty(boxwidth);
            boxwidth=boxwdefault;
        end
    else
        if isempty(boxwidth);
            boxwidth=boxwdefault;
%         else
%             boxwdefault=boxwidth;
        end
        if a1==1;
            boxwidth=abs(boxwdefault*(xvals(2)-xvals(1)));
        elseif a1==ngroups;
            boxwidth=abs(boxwdefault*(xvals(a1)-xvals(a1-1)));
        else
            boxwidth=boxwdefault*min([abs(xvals(a1+1)-xvals(a1)),abs(xvals(a1)-xvals(a1-1))]);
        end
    end
    xvec=[xvals(a1)-boxwidth/2 xvals(a1)+boxwidth/2 xvals(a1)+boxwidth/2 xvals(a1)-boxwidth/2];
    yvec=[pctles(1) pctles(1) pctles(2) pctles(2)];
    if size(colormatrix,1)==ngroups;
        colorvec=colormatrix(a1,:);
    else
        colorvec=colormatrix(1,:);
    end
    if size(boxfillcolormat,1)==ngroups;
        boxfillcolor=boxfillcolormat(a1,:);
    else
        boxfillcolor=boxfillcolormat(1,:);
    end
    if sdlength>0;
        DoErrorBarsCvec(xvals(a1),m,(m+sd),0.3*boxwidth,linw1,colorvec);
        DoErrorBarsCvec(xvals(a1),m,(m-sd),0.3*boxwidth,linw1,colorvec);
        hold on;
    end
%     hh=patch(xvec,yvec,'edgecolor',colorvec,'faccecolor',boxfillcolor,'linewidth',linw1);

    hh=patch(xvec,yvec,boxfillcolor,'edgecolor',colorvec,'linewidth',linw1);
    if nofacecolor
        set(hh,'facecolor','none');
    end
    hold on;
    line(xvec(1:2),[med med],'color',colorvec,'linewidth',linw1);
    if meanmarksize>0;
        hold on;
        plot(xvals(a1),m,'o','markeredgecolor',colorvec,'markerfacecolor',colorvec,'linewidth',linw,'markersize',meanmarksize);
    end
end
