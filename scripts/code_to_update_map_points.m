 
 
 %%
 
 left_labels = {'AF7','F5','F3','FC3','C3','CP3','P3','P5','PO7','PO9'};
 right_labels = {'AF8','F6','F4','FC4','C4','CP4','P4','P6','PO8','PO10'};
 top_left_labels = {'FP1','AF3','F1','FC1','C1','CP1','P1','PO3','O1','I1'};
 top_right_labels = {'FP2','AF4','F2','FC2','C2','CP2','P2','PO4','O2','I2'};
 
 load('latticeExtensions.mat');
 
 %%
 load('map-sideLeft.mat')
 figure; plot(outline(:, 1), outline(:, 2)); hold on; 
%  plot(points(:, 1), points(:, 2),'*');
 for u= 1:size(points,1)
     text(points(u, 1), points(u, 2),labels{u});
 end
 hold off


eSeamXpts = latticeExtensions.lefteSeamXpts(:,1:2);
ii = find(eSeamXpts(:,1) == 0 & eSeamXpts(:,2) == 0);
SeamXpts = eSeamXpts;
SeamXpts(ii,:) = [];

idx = find(ismember(labels, left_labels )==1);
load('state.mat','sideLeftOutline');

pts_add =[];
labels_add = labels(idx);
for u=1:length(idx)
    pt = points(idx(u),:)-[min(sideLeftOutline(:,1)) min(sideLeftOutline(:,2))];
    dist = sqrt(sum((SeamXpts-pt).^2,2));
    [~, n_idx] = min(dist);
    n_idx = find(eSeamXpts(:,1) == SeamXpts(n_idx,1) & eSeamXpts(:,2) == SeamXpts(n_idx,2));
    Movingpts = [latticeExtensions.lefteSeamXpts(n_idx,1) latticeExtensions.lefteSeamXpts(n_idx,2); latticeExtensions.lefteSeamXpts(n_idx,5) latticeExtensions.lefteSeamXpts(n_idx,6)];
    FixedPts = [latticeExtensions.toplefteSeamXpts(n_idx,1) latticeExtensions.toplefteSeamXpts(n_idx,2); latticeExtensions.toplefteSeamXpts(n_idx,5) latticeExtensions.toplefteSeamXpts(n_idx,6)];
    tform = fitgeotrans(Movingpts,FixedPts,'nonreflectivesimilarity');
    [Txpt,Typt] = transformPointsForward(tform, pt(1), pt(2));
    pts_add = [pts_add; [Txpt Typt]];
end

load('map-top.mat')
points = [points; pts_add];
for u=1:length(labels_add)
    labels{end+1} = labels_add{u};
end
save('updated_map-top.mat','labels','points','outline','panel');
 figure; plot(outline(:, 1), outline(:, 2)); hold on; 
%  plot(points(:, 1), points(:, 2),'*');
 for u= 1:size(points,1)
     text(points(u, 1), points(u, 2),labels{u});
 end
 for u=1:size(pts_add)
     text(pts_add(u, 1), pts_add(u, 2),labels_add{u});
 end
 hold off
 axis image
%%

 load('map-top.mat')
 figure; plot(outline(:, 1), outline(:, 2)); hold on; 
%  plot(points(:, 1), points(:, 2),'*');
 for u= 1:size(points,1)
     text(points(u, 1), points(u, 2),labels{u});
 end
 hold off


eSeamXpts = latticeExtensions.toplefteSeamXpts(:,1:2);
ii = find(eSeamXpts(:,1) == 0 & eSeamXpts(:,2) == 0);
SeamXpts = eSeamXpts;
SeamXpts(ii,:) = [];

idx = find(ismember(labels, top_left_labels )==1);
load('state.mat','topOutline');

pts_add =[];
labels_add = labels(idx);
for u=1:length(idx)
    pt = points(idx(u),:)-[min(topOutline(:,1)) min(topOutline(:,2))];
    dist = sqrt(sum((SeamXpts-pt).^2,2));
    [~, n_idx] = min(dist);
    n_idx = find(eSeamXpts(:,1) == SeamXpts(n_idx,1) & eSeamXpts(:,2) == SeamXpts(n_idx,2));
    Movingpts = [latticeExtensions.toplefteSeamXpts(n_idx,1) latticeExtensions.toplefteSeamXpts(n_idx,2); latticeExtensions.toplefteSeamXpts(n_idx,5) latticeExtensions.toplefteSeamXpts(n_idx,6)];
    FixedPts = [latticeExtensions.lefteSeamXpts(n_idx,1) latticeExtensions.lefteSeamXpts(n_idx,2); latticeExtensions.lefteSeamXpts(n_idx,5) latticeExtensions.lefteSeamXpts(n_idx,6)];
    tform = fitgeotrans(Movingpts,FixedPts,'nonreflectivesimilarity');
    [Txpt,Typt] = transformPointsForward(tform, pt(1), pt(2)); 
    pts_add = [pts_add; [Txpt Typt]+[min(sideLeftOutline(:,1)) min(sideLeftOutline(:,2))];];
%     labels_add{u} = labels{idx}
end


load('map-sideLeft.mat')
points = [points; pts_add];
for u=1:length(labels_add)
    labels{end+1} = labels_add{u};
end
save('updated_map-sideLeft.mat','labels','points','outline','panel');
 figure; plot(outline(:, 1), outline(:, 2)); hold on; 
%  plot(points(:, 1), points(:, 2),'*');
 for u= 1:size(points,1)
     text(points(u, 1), points(u, 2),labels{u});
 end
 for u=1:size(pts_add)
     text(pts_add(u, 1), pts_add(u, 2),labels_add{u});
 end
 hold off
% axis image
%%
 load('map-sideRight.mat')
 figure; plot(outline(:, 1), outline(:, 2)); hold on; 
%  plot(points(:, 1), points(:, 2),'*');
 for u= 1:size(points,1)
     text(points(u, 1), points(u, 2),labels{u});
 end
 hold off


eSeamXpts = latticeExtensions.righteSeamXpts(:,1:2);
ii = find(eSeamXpts(:,1) == 0 & eSeamXpts(:,2) == 0);
SeamXpts = eSeamXpts;
SeamXpts(ii,:) = [];

idx = find(ismember(labels, right_labels )==1);
load('state.mat','sideRightOutline');

pts_add =[];
labels_add = labels(idx);
for u=1:length(idx)
    pt = points(idx(u),:)-[min(sideRightOutline(:,1)) min(sideRightOutline(:,2))];
    dist = sqrt(sum((SeamXpts-pt).^2,2));
    [~, n_idx] = min(dist);
    n_idx = find(eSeamXpts(:,1) == SeamXpts(n_idx,1) & eSeamXpts(:,2) == SeamXpts(n_idx,2));
    Movingpts = [latticeExtensions.righteSeamXpts(n_idx,1) latticeExtensions.righteSeamXpts(n_idx,2); latticeExtensions.righteSeamXpts(n_idx,5) latticeExtensions.righteSeamXpts(n_idx,6)];
    FixedPts = [latticeExtensions.toprighteSeamXpts(n_idx,1) latticeExtensions.toprighteSeamXpts(n_idx,2); latticeExtensions.toprighteSeamXpts(n_idx,5) latticeExtensions.toprighteSeamXpts(n_idx,6)];
    tform = fitgeotrans(Movingpts,FixedPts,'nonreflectivesimilarity');
    [Txpt,Typt] = transformPointsForward(tform, pt(1), pt(2));
    pts_add = [pts_add; [Txpt Typt]];
%     labels_add{u} = labels{idx}
end


load('updated_map-top.mat')
points = [points; pts_add];
for u=1:length(labels_add)
    labels{end+1} = labels_add{u};
end
save('updated_map-top.mat','labels','points','outline','panel');
 figure; plot(outline(:, 1), outline(:, 2)); hold on; 
%  plot(points(:, 1), points(:, 2),'*');
 for u= 1:size(points,1)
     text(points(u, 1), points(u, 2),labels{u});
 end
 for u=1:size(pts_add)
     text(pts_add(u, 1), pts_add(u, 2),labels_add{u});
 end
 hold off
 axis image
 
 %%
 
 load('map-top.mat')
 figure; plot(outline(:, 1), outline(:, 2)); hold on; 
%  plot(points(:, 1), points(:, 2),'*');
 for u= 1:size(points,1)
     text(points(u, 1), points(u, 2),labels{u});
 end
 hold off


eSeamXpts = latticeExtensions.toprighteSeamXpts(:,1:2);
ii = find(eSeamXpts(:,1) == 0 & eSeamXpts(:,2) == 0);
SeamXpts = eSeamXpts;
SeamXpts(ii,:) = [];

idx = find(ismember(labels, top_right_labels )==1);
load('state.mat','topOutline');

pts_add =[];
labels_add = labels(idx);
for u=1:length(idx)
    pt = points(idx(u),:)-[min(topOutline(:,1)) min(topOutline(:,2))];
    dist = sqrt(sum((SeamXpts-pt).^2,2));
    [~, n_idx] = min(dist);
    n_idx = find(eSeamXpts(:,1) == SeamXpts(n_idx,1) & eSeamXpts(:,2) == SeamXpts(n_idx,2));
    Movingpts = [latticeExtensions.toprighteSeamXpts(n_idx,1) latticeExtensions.toprighteSeamXpts(n_idx,2); latticeExtensions.toprighteSeamXpts(n_idx,5) latticeExtensions.toprighteSeamXpts(n_idx,6)];
    FixedPts = [latticeExtensions.righteSeamXpts(n_idx,1) latticeExtensions.righteSeamXpts(n_idx,2); latticeExtensions.righteSeamXpts(n_idx,5) latticeExtensions.righteSeamXpts(n_idx,6)];
    tform = fitgeotrans(Movingpts,FixedPts,'nonreflectivesimilarity');
    [Txpt,Typt] = transformPointsForward(tform, pt(1), pt(2)); 
    pts_add = [pts_add; [Txpt Typt]+[min(sideRightOutline(:,1)) min(sideRightOutline(:,2))];];
%     labels_add{u} = labels{idx}
end


load('map-sideRight.mat')
points = [points; pts_add];
for u=1:length(labels_add)
    labels{end+1} = labels_add{u};
end
save('updated_map-sideRight.mat','labels','points','outline','panel');
 figure; plot(outline(:, 1), outline(:, 2)); hold on; 
%  plot(points(:, 1), points(:, 2),'*');
 for u= 1:size(points,1)
     text(points(u, 1), points(u, 2),labels{u});
 end
 for u=1:size(pts_add)
     text(pts_add(u, 1), pts_add(u, 2),labels_add{u});
 end
 hold off