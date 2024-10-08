function [STLcoords] = getStlCoordinates(grommets, holders, aux, sideOutlineL, sideOutlineR, sideLeftPanel, topOutline, NsideOutsidePoints, IDXchStrapPoints)
%% gives back absolute coordinates for STL files to enable blender to paste elements at the right positions

%% mirror sideOutline points to get sideLeft coordinates, based on what reflectX does for poly
% find centroid using the poly function, then mirror the sideOutline
% [x, y] = boundingbox(sideLeftPanel);

% x = [min(sideOutline(:,1)), max(sideOutline(:,1))];
% y = [min(sideOutline(:,2)), max(sideOutline(:,2))];

% cS = [mean(x) mean(y)];
% coords.topGrommets = cat(1, grommets(arrayfun(@(x) strcmp(x.panel, 'top') && ~strcmp(x.type, '#DUMMY') && ~ismember('short-separation', x.flags), grommets)).posPanel);

%% Grommets
% right side
% STLcoords.sideRightPos = cat(1, grommets(arrayfun(@(x) strcmp(x.panel, 'sideRight') && ~strcmp(x.type, '#DUMMY') && ~ismember('short-separation', x.flags), grommets)).posPanel);
% STLcoords.sideRightID  = cat(1, grommets(arrayfun(@(x) strcmp(x.panel, 'sideRight') && ~strcmp(x.type, '#DUMMY') && ~ismember('short-separation', x.flags), grommets)).type);
% STLcoords.sideRightRot = cat(1, grommets(arrayfun(@(x) strcmp(x.panel, 'sideRight') && ~strcmp(x.type, '#DUMMY') && ~ismember('short-separation', x.flags),  grommets)).rot);
% % left side
% STLcoords.sideLeftPos = cat(1, grommets(arrayfun(@(x) strcmp(x.panel, 'sideLeft') && ~strcmp(x.type, '#DUMMY') && ~ismember('short-separation', x.flags), grommets)).posPanel);
% STLcoords.sideLeftID  = cat(1, grommets(arrayfun(@(x) strcmp(x.panel, 'sideLeft') && ~strcmp(x.type, '#DUMMY') && ~ismember('short-separation', x.flags), grommets)).type);
% STLcoords.sideLeftRot  = cat(1, grommets(arrayfun(@(x) strcmp(x.panel, 'sideLeft') && ~strcmp(x.type, '#DUMMY') && ~ismember('short-separation', x.flags), grommets)).rot);
% % top pieces
% STLcoords.top1Pos = cat(1, grommets(arrayfun(@(x) strcmp(x.panel, 'top') && x.panelIndex ==1 && ~strcmp(x.type, '#DUMMY') && ~ismember('short-separation', x.flags), grommets)).posPanel);
% STLcoords.top1ID  = cat(1, grommets(arrayfun(@(x) strcmp(x.panel, 'top') && ~strcmp(x.type, '#DUMMY') && x.panelIndex ==1 && ~ismember('short-separation', x.flags),  grommets)).type);
% STLcoords.top1Rot  = cat(1, grommets(arrayfun(@(x) strcmp(x.panel, 'top') && ~strcmp(x.type, '#DUMMY') && x.panelIndex ==1 && ~ismember('short-separation', x.flags), grommets)).rot);
% STLcoords.top2Pos = cat(1, grommets(arrayfun(@(x) strcmp(x.panel, 'top') && x.panelIndex ==2 && ~strcmp(x.type, '#DUMMY') && ~ismember('short-separation', x.flags), grommets)).posPanel);
% STLcoords.top2ID  = cat(1, grommets(arrayfun(@(x) strcmp(x.panel, 'top') && ~strcmp(x.type, '#DUMMY') && x.panelIndex ==2 && ~ismember('short-separation', x.flags), grommets)).type);
% STLcoords.top2Rot  = cat(1, grommets(arrayfun(@(x) strcmp(x.panel, 'top') && ~strcmp(x.type, '#DUMMY') && x.panelIndex ==2 && ~ismember('short-separation', x.flags), grommets)).rot);

%% Grommets
% right side
% STLcoords.sideRightPos = cat(1, grommets(arrayfun(@(x) strcmp(x.panel, 'sideRight') && ~strcmp(x.type, '#DUMMY') && ~ismember('short-separation', x.flags), grommets)).posPanel);
STLcoords.sideRightPos = cat(1, grommets(arrayfun(@(x) strcmp(x.panel, 'sideRight') && ~strcmp(x.type, '#DUMMY'), grommets)).posPanel);
STLcoords.sideRightID  = {grommets(arrayfun(@(x) strcmp(x.panel, 'sideRight') && ~strcmp(x.type, '#DUMMY'), grommets)).type};
STLcoords.sideRightRot = cat(1, grommets(arrayfun(@(x) strcmp(x.panel, 'sideRight') && ~strcmp(x.type, '#DUMMY'),  grommets)).rot);
% left side
% % STLcoords.sideLeftPos = cat(1, grommets(arrayfun(@(x) strcmp(x.panel, 'sideLeft') && ~strcmp(x.type, '#DUMMY') && ~ismember('short-separation', x.flags), grommets)).posPanel);
STLcoords.sideLeftPos = cat(1, grommets(arrayfun(@(x) strcmp(x.panel, 'sideLeft') && ~strcmp(x.type, '#DUMMY'), grommets)).posPanel);
STLcoords.sideLeftID  = {grommets(arrayfun(@(x) strcmp(x.panel, 'sideLeft') && ~strcmp(x.type, '#DUMMY'), grommets)).type};
STLcoords.sideLeftRot  = cat(1, grommets(arrayfun(@(x) strcmp(x.panel, 'sideLeft') && ~strcmp(x.type, '#DUMMY'), grommets)).rot);
% top pieces
% STLcoords.top1Pos = cat(1, grommets(arrayfun(@(x) strcmp(x.panel, 'top') && x.panelIndex ==1 && ~strcmp(x.type, '#DUMMY') && ~ismember('short-separation', x.flags), grommets)).posPanel);
STLcoords.top1Pos = cat(1, grommets(arrayfun(@(x) strcmp(x.panel, 'top') && x.panelIndex ==1 && ~strcmp(x.type, '#DUMMY'), grommets)).posPanel);
STLcoords.top1ID  = {grommets(arrayfun(@(x) strcmp(x.panel, 'top') && ~strcmp(x.type, '#DUMMY') && x.panelIndex ==1, grommets)).type};
STLcoords.top1Rot  = cat(1, grommets(arrayfun(@(x) strcmp(x.panel, 'top') && ~strcmp(x.type, '#DUMMY') && x.panelIndex ==1, grommets)).rot);
% STLcoords.top2Pos = cat(1, grommets(arrayfun(@(x) strcmp(x.panel, 'top') && x.panelIndex ==2 && ~strcmp(x.type, '#DUMMY') && ~ismember('short-separation', x.flags), grommets)).posPanel);
STLcoords.top2Pos = cat(1, grommets(arrayfun(@(x) strcmp(x.panel, 'top') && x.panelIndex ==2 && ~strcmp(x.type, '#DUMMY'), grommets)).posPanel);
% STLcoords.top2ID  = cat(1, grommets(arrayfun(@(x) strcmp(x.panel, 'top') && ~strcmp(x.type, '#DUMMY') && x.panelIndex ==2, grommets)).type);
STLcoords.top2ID  = {grommets(arrayfun(@(x) strcmp(x.panel, 'top') && ~strcmp(x.type, '#DUMMY') && x.panelIndex ==2, grommets)).type};
STLcoords.top2Rot  = cat(1, grommets(arrayfun(@(x) strcmp(x.panel, 'top') && ~strcmp(x.type, '#DUMMY') && x.panelIndex ==2, grommets)).rot);



%% Earpoints
% calculate midpoint for ear using first and last point on lower
% sideOutline from SVG file and an empirical factor sfact along the line
% sfact = 3/5;
% earrefR = sideOutlineR([end-NsideOutsidePoints,end],:);
% STLcoords.sideRightEar =  sfact*earrefR(2,:)+(1-sfact)*earrefR(1,:);
% 
% % do L as well
% earrefL = sideOutlineL([end-NsideOutsidePoints,end],:);
% STLcoords.sideLeftEar =  sfact*earrefL(2,:)+(1-sfact)*earrefL(1,:);

%% ChinStrap
% set chin strap holders at midpoint of given coordinates
STLcoords.sideRightChStrap =  (sideOutlineR(end-IDXchStrapPoints(1),:)+sideOutlineR(end-IDXchStrapPoints(2),:))/2;
STLcoords.sideLeftChStrap =  (sideOutlineL(end-IDXchStrapPoints(1),:)+sideOutlineL(end-IDXchStrapPoints(2),:))/2;

STLcoords.sideRightChStrap_2 =  sideOutlineR(end-IDXchStrapPoints(1)+3,:);
STLcoords.sideLeftChStrap_2 =  sideOutlineL(end-IDXchStrapPoints(1)+3,:);





