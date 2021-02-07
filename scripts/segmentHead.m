function [sideLeftOutline, sideRightOutline, outline, vertSideLeft, vertSideRight, vertTop] = segmentHead(vHead, fHead, refpts, side_piece)
%SEGMENTHEAD
%   This function was designed to take in a scaled head model and create a
%   set of outlines that will be used to create the cap. The function finds
%   two thresholds that will be used as the borders of the outlines. The
%   different panels of the cap will be joined along these threshold lines.

%% parameters
center_size = 3 / 7; % width of center segment
extend_back = 0; % mm, how far below (or above) inion
extend_front = -12.5; % mm, how far below (or more likely, above) nasion

%% visualize head
figure(1);
subplot(1, 3, 1);
h = trisurf(fHead, vHead(:, 1), vHead(:, 2), vHead(:, 3));
set(h, 'LineStyle', 'none');
set(h, 'FaceAlpha', 0.5);
xlabel('x'); ylabel('y'); zlabel('z');
axis image;

%% get z-axis
% midline plane: use Nz, Cz and Iz
Nz = getReferencePt(refpts, 'Nz');
Cz = getReferencePt(refpts, 'Cz');
Iz = getReferencePt(refpts, 'Iz');
plane_pts_midline = [Nz; Cz; Iz];
plane_midline = makePlane(plane_pts_midline);

slice_midline = slice3D(vHead, fHead, plane_midline);
[~, pts] = projectPlane(plane_midline, plane_pts_midline);

%% segment head
% segment the head surface into three pieces

% get left right points
edge_right = getReferencePt(refpts, 'RPA');
edge_left = getReferencePt(refpts, 'LPA');

% sometimes backwards...
if planeDistance(plane_midline, edge_right) < 0
    [edge_left, edge_right] = deal(edge_right, edge_left);
end

% find width of center section
center_width = norm(edge_left - edge_right) * center_size;

% find norm
norm_x = planeNormal(plane_midline);

% get sign for right edge and make right edge plane
s = sign(dot(norm_x, edge_right - plane_midline(1:3)));
plane_right = plane_midline + [(s * norm_x * center_width / 2) 0 0 0 0 0 0];

% get sign for left edge and make left edge plane
s = sign(dot(norm_x, edge_left - plane_midline(1:3)));
plane_left = plane_midline + [(s * norm_x * center_width / 2) 0 0 0 0 0 0];

% get slices
slice_right = slice3D(vHead, fHead, plane_right);
slice_left = slice3D(vHead, fHead, plane_left);

% move Nz and Iz to form the front and back of the cap (in 3D coordinates)
front = projectFromPlane(plane_midline, sliceMovePoint(slice_midline, pts(1, :), pts(2, :), extend_front));
back = projectFromPlane(plane_midline, sliceMovePoint(slice_midline, pts(3, :), pts(2, :), extend_back));

% get path
[~, p] = projectPlane(plane_right, [front; Cz; back]);
edge_right = sliceGetEdge(slice_right, p(1, :), p(3, :), p(2, :));
[~, p] = projectPlane(plane_left, [front; Cz; back]);
edge_left = sliceGetEdge(slice_left, p(1, :), p(3, :), p(2, :));

%% get left width length
length_left = sum(sqrt(sum((edge_left(2:end, :) - edge_left(1:(end - 1), :)) .^ 2, 2)));

%% get top width
length_width = segmentTopWidth(vHead, fHead, front, plane_left, edge_left, plane_right, edge_right);

%% make top outline

% outline
outline = [0 0; ... % back left
    0 length_left; ... % front left
    length_width length_left; ... % front right
    length_width 0]; % back right

% vertices
vertTop = [projectFromPlane(plane_left, [edge_left(end, :); edge_left(1, :)]); ...
    projectFromPlane(plane_right, [edge_right(1, :); edge_right(end, :)])];

figure(1);
subplot(1, 3, 1);
hold on;
scatter3(vertTop(:, 1), vertTop(:, 2), vertTop(:, 3));
hold off;

subplot(1, 3, 2);
scatter(outline(:, 1), outline(:, 2));
axis image;

%% make LEFT side outline

% get points for transform
points_fixed = [edge_left(1, :); edge_left(end, :)];
points_moving = [side_piece(1, :); side_piece(end, :)];

% make transform
tform_cpp = fitgeotrans(points_moving, points_fixed, 'nonreflectivesimilarity');

% apply transform
tmp = transformPointsForward(tform_cpp, side_piece);
tmp = flipud(tmp);

% vertices
vertSideLeft = [nan(size(tmp, 1), size(tmp, 2) + 1); projectFromPlane(plane_left, edge_left)];

% outline
sideLeftOutline = [tmp; edge_left];

% reflect
sideLeftOutline = reflectXOutline(sideLeftOutline);

%% make RIGHT side outline

% get points for transform
points_fixed = [edge_right(1, :); edge_right(end, :)];
points_moving = [side_piece(1, :); side_piece(end, :)];

% make transform
tform_cpp = fitgeotrans(points_moving, points_fixed, 'nonreflectivesimilarity');

% apply transform
tmp = transformPointsForward(tform_cpp, side_piece);
tmp = flipud(tmp);

% vertices
vertSideRight = [nan(size(tmp, 1), size(tmp, 2) + 1); projectFromPlane(plane_right, edge_right)];

% outline
sideRightOutline = [tmp; edge_right];

end