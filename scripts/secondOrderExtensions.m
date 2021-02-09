function [includeLeftExtension, includeRightExtension, includeTopLeftExtension, includeTopRightExtension] = secondOrderExtensions(grommets, C)
leftSideIdx = [];
rightSideIdx = [];
topIdx = [];
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
leftPos = [grommets(leftSideIdx).posPanel];
rightPos = [grommets(rightSideIdx).posPanel];
topPos = [grommets(topIdx).posPanel];
dist_lessthan = 25/2;
%%
includeLeftExtension = zeros(size(C(2).eSeamExt,1),1);
for u = 1:2:length(leftPos)
    pt = leftPos(u:u+1);
    dist = sqrt(sum((C(2).eSeamExt(:,1:2)-pt).^2,2));
    includeLeftExtension(dist<=dist_lessthan) = 1;
    pt = leftPos(u:u+1);
    dist = sqrt(sum((C(2).eSeamExt(:,3:4)-pt).^2,2));
    includeLeftExtension(dist<=dist_lessthan) = 1;
end
%%
includeRightExtension = zeros(size(C(3).eSeamExt,1),1);
for u = 1:2:length(rightPos)
    pt = rightPos(u:u+1);
    dist = sqrt(sum((C(3).eSeamExt(:,1:2)-pt).^2,2));
    includeRightExtension(dist<=dist_lessthan) = 1;
    pt = rightPos(u:u+1);
    dist = sqrt(sum((C(3).eSeamExt(:,3:4)-pt).^2,2));
    includeRightExtension(dist<=dist_lessthan) = 1;
end
%%
includeTopLeftExtension = zeros(size(C(2).eSeamExtRot,1),1);
for u = 1:2:length(topPos)
    pt = topPos(u:u+1);
    dist = sqrt(sum((C(2).eSeamExtRot(:,1:2)-pt).^2,2));
    includeTopLeftExtension(dist<=dist_lessthan) = 1;
    pt = topPos(u:u+1);
    dist = sqrt(sum((C(2).eSeamExtRot(:,3:4)-pt).^2,2));
    includeTopLeftExtension(dist<=dist_lessthan) = 1;
end
%%
includeTopRightExtension = zeros(size(C(3).eSeamExtRot,1),1);
for u = 1:2:length(topPos)
    pt = topPos(u:u+1);
    dist = sqrt(sum((C(3).eSeamExtRot(:,1:2)-pt).^2,2));
    includeTopRightExtension(dist<=dist_lessthan) = 1;
    pt = topPos(u:u+1);
    dist = sqrt(sum((C(3).eSeamExtRot(:,3:4)-pt).^2,2));
    includeTopRightExtension(dist<=dist_lessthan) = 1;
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