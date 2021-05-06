function [includeLeftExtension, includeRightExtension, includeTopLeftExtension, includeTopRightExtension] = secondOrderExtensions(grommets, C)
leftSideIdx = [];
rightSideIdx = [];
topIdx = [];
load('grommetDimensions.mat','grommetDimensions')
for u = 1:length(grommets)
    if strcmp(grommets(u).panel,'sideLeft')
        leftSideIdx = [leftSideIdx; u];
    elseif strcmp(grommets(u).panel,'sideRight')
        rightSideIdx = [rightSideIdx; u];
    elseif strcmp(grommets(u).panel,'top')
        topIdx = [topIdx; u];
    end
end
%%
includeLeftExtension = zeros(size(C(2).eSeamExt,1),1);
for u = 1:length(leftSideIdx)
    pt = grommets(leftSideIdx(u)).posPanel;
    grommet_type = grommets(leftSideIdx(u)).type;
    if ~strcmp(grommet_type,'#DUMMY')
        idx = find(ismember(grommetDimensions(:,1),grommet_type)==1);
        grommet_max_size = max(grommetDimensions{idx,2:3})/2;
        dist = sqrt(sum((C(2).eSeamExt(:,1:2)-pt).^2,2));
        includeLeftExtension(dist<=grommet_max_size) = 1;
        dist = sqrt(sum((C(2).eSeamExt(:,3:4)-pt).^2,2));
        includeLeftExtension(dist<=grommet_max_size) = 1;
    end
end
%%
includeRightExtension = zeros(size(C(3).eSeamExt,1),1);
for u = 1:length(rightSideIdx)
    pt = grommets(rightSideIdx(u)).posPanel;
    grommet_type = grommets(rightSideIdx(u)).type;
    if ~strcmp(grommet_type,'#DUMMY')
        idx = find(ismember(grommetDimensions(:,1),grommet_type)==1);
        grommet_max_size = max(grommetDimensions{idx,2:3})/2;
        dist = sqrt(sum((C(3).eSeamExt(:,1:2)-pt).^2,2));
        includeRightExtension(dist<=grommet_max_size) = 1;
        dist = sqrt(sum((C(3).eSeamExt(:,3:4)-pt).^2,2));
        includeRightExtension(dist<=grommet_max_size) = 1;
    end
end
%%
includeTopLeftExtension = zeros(size(C(2).eSeamExtRot,1),1);
includeTopRightExtension = zeros(size(C(3).eSeamExtRot,1),1);
for u = 1:length(topIdx)
    pt = grommets(topIdx(u)).posPanel;
    grommet_type = grommets(topIdx(u)).type;
    if ~strcmp(grommet_type,'#DUMMY')
        idx = find(ismember(grommetDimensions(:,1),grommet_type)==1);
        grommet_max_size = max(grommetDimensions{idx,2:3})/2;
        dist = sqrt(sum((C(2).eSeamExtRot(:,1:2)-pt).^2,2));
        includeTopLeftExtension(dist<=grommet_max_size) = 1;
        dist = sqrt(sum((C(3).eSeamExtRot(:,1:2)-pt).^2,2));
        includeTopRightExtension(dist<=grommet_max_size) = 1;
        dist = sqrt(sum((C(2).eSeamExtRot(:,3:4)-pt).^2,2));
        includeTopLeftExtension(dist<=grommet_max_size) = 1;
        dist = sqrt(sum((C(3).eSeamExtRot(:,3:4)-pt).^2,2));
        includeTopRightExtension(dist<=grommet_max_size) = 1;
    end
end

%%
% if length of the overlap strut is smaller than add second order extension
for u = 1:size(C(2).eSeamExt,1)
    ext_length = sqrt(sum((C(2).eSeamExt(u,1:2)-C(2).eSeamExt(u,3:4)).^2));
    if ext_length <= 7
        includeLeftExtension(u) = 1;
    end
end


for u = 1:size(C(3).eSeamExt,1)
    ext_length = sqrt(sum((C(3).eSeamExt(u,1:2)-C(3).eSeamExt(u,3:4)).^2));
    if ext_length <= 7
        includeRightExtension(u) = 1;
    end
end
