% modify the map-panel-10_5.mat files to have only refpts INSIDE each panel
% 
% take the template/state.mat panelOutline and panel3DVertices and make
% uniform spacing of 1mm between vertices
% Do this within the code... no need to update state.mat
% Note that this solves the problem of only 4 vertices we saw in the top
% panel.
%
% modify mapGrommets to include distances to every Nth vertices amongst
% distances to refpts. Consider the closest 6. Nth can be 10 or maybe 7.
%         distancesOutline = ones(num2, 1)*1e9;
%         for i = 1:6:num2
%             distancesOutline(i) = surfaceDist(outline3DVertices(i, :), grommet);
%         end



foo=load('template/state.mat');

%% Load side piece
% New function reads in outline from a SVG file
fname = 'svg/side_piece2.svg';
pathid = 'side_piece';
% side panel outline and number of side panel points that define the outside vertices
[side_piece, sideOutsidePoints] = svgImport(fname, pathid);


[sideLeftOutline, sideRightOutline, topOutline, sideLeft3DVertices, sideRight3DVertices, top3DVertices] = segmentHead(foo.vHead, foo.fHead, foo.refpts, side_piece);

% SK - MAAKE UNIFORM SPACING OF VERTICES... 1 MM

figure(2)
plot(sideLeftOutline(:,1), sideLeftOutline(:,2) ); 
axis image

hold on
plot(foo.sideLeftOutline(:,1), foo.sideLeftOutline(:,2), 'r.' ); 
hold off