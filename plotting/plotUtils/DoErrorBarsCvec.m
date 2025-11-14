function DoErrorBarsCvec(x,ystart,yend,width,linewidth,linecol);

if isempty(x);
    x=1:length(ystart);
end
if isempty(width);
    width=0.3*mean(diff(x));
end
if isempty(linewidth);
    linewidth=1;
end
if isempty(linecol);
    linecol=[0 0 0];
end

if (size(ystart,1)>size(ystart,2));
    ystart=ystart';
end
if (size(yend,1)>size(yend,2));
    yend=yend';
end
if (size(x,1)>size(x,2));
    x=x';
end

line([x;x],[ystart;yend],'linewidth',linewidth,'color',linecol);
line([(x-width/2);(x+width/2)],[yend;yend],'linewidth',linewidth,'color',linecol);

