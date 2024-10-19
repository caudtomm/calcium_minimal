function anat_regions = selectAnatRegions(subject, ignore_previous)
arguments
    subject Subject
    ignore_previous logical = false
end
%%% REGION SELECTION

anatomy = subject.reference_img;


%%

a = subject.ROImap;
rois = sort(unique(a));
centroids = [];
for i_roi = 2:numel(unique(rois))
    c = double(a==rois(i_roi));
    b = regionprops(c,'Centroid');
    d = cat(1,b.Centroid);
    if size(d,1)>1; dbstop; end
    centroids = [centroids; d];
end; clear a b c d i_roi


figure(6384790)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

a = subject.ROImap; a(a>0) = quantile(anatomy(:), .1);
imagesc(anatomy + a); colormap('gray'); clear a
hold on; scatter(centroids(:,1),centroids(:,2),'g')


% region selection

%initialize vars
anat_regions = subject.anat_regions;
c = optimalcolors(11); c(1,:)=[];
if ignore_previous || isempty(subject.anat_regions)
    i_region = 1;
    isInROI = [];
    [h, cells_in_region, region_names] = deal({});
else
    h = anat_regions.position;
    cells_in_region = anat_regions.cells(1:end-1);
    region_names = anat_regions.names(1:end-1);
    
    isInROI = zeros(size(centroids,1),1);
    for i_cell = 1:numel(isInROI)
        for i_region = 1:numel(region_names)
            itis = ismember(rois(i_cell+1),cells_in_region{i_region});
            if itis
                isInROI(i_cell) = 1;
                continue
            end
        end
    end
    
    %show preexisting regions
    for i_region = 1:numel(region_names)
        h1 = drawpolygon('FaceAlpha',0.1, ...
                         'Color',c(i_region,:), ...
                         'InteractionsAllowed','none', ...
                         'Label',region_names{i_region}, ...
                         'LabelAlpha',0, ...
                         'LabelTextColor',c(i_region,:), ...
                         'Position', h{i_region}, ...
                         'Deletable', false);
    end
    i_region = i_region+1;
end

%execute prompt
while true
    prompt = ['Name of region #',num2str(i_region),': '];
    region_name = input(prompt,'s');
    
    while true
        disp('draw!')
        
        h1 = drawpolygon('FaceAlpha',0.1, ...
                         'Color',c(i_region,:), ...
                         'InteractionsAllowed','reshape', ...
                         'Label',region_name, ...
                         'LabelAlpha',0, ...
                         'LabelTextColor',c(i_region,:));
        
        is_keep = true;
        while true
            s = input('keep region? ','s');

            if strcmp(s,'y')||strcmp(s,'Y')
                break
            elseif strcmp(s,'n')||strcmp(s,'N')
                is_keep = false;
                break               
            else
                disp('Input invalid')
            end
        end
        
        if is_keep
            h{end+1} = h1.Position;
            break
        else
            delete(h1)
        end
    end
    
    tf = inROI(h1,centroids(:,1),centroids(:,2));
    isInROI = [isInROI, tf];
    idx = find(tf);
    cells_in_region{end+1} = rois(idx+1);
    region_names{end+1} = region_name;
    
    i_region = i_region+1;
    
    is_goon = true;
    while true
        s = input('done? ','s');

        if strcmp(s,'y')||strcmp(s,'Y')
            is_goon = false;
            break
        elseif strcmp(s,'n')||strcmp(s,'N')
            break
        else
            disp('Input invalid')
        end
    end
    
    if ~is_goon; break; end
end

% add 'other' region
a = sum(isInROI,2);
idx = find(a==0);
region_names{end+1} = 'other';
cells_in_region{end+1} = rois(idx+1);


%% checking consistency

% isInROI
a = sum(isInROI,2);
b = find(a>1);
if ~isempty(b)
    disp([num2str(numel(b)),' cells show multiple assignments!'])
    disp('this program will do nothing about it! (:')
end


%% saving to data struct


anat_regions.position = h;
anat_regions.names = region_names;
anat_regions.cells = cells_in_region;


close(6384790)






















