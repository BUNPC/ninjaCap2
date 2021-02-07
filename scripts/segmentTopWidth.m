function length_width = segmentTopWidth(vHead, fHead, front, plane_left, edge_left, plane_right, edge_right)
%SEGMENTTOPWIDTH Summary of this function goes here
%   Detailed explanation goes here

% make points
front_points = [...
    projectFromPlane(plane_left, edge_left(1, :)); ...
    front; ...
    projectFromPlane(plane_right, edge_right(1, :)) ...
    ];

% make plane
plane = makePlane(front_points);

% get path around
slice = slice3D(vHead, fHead, plane);

% project points
[~, front_points_2d] = projectPlane(plane, front_points);

% get length
length_width = sliceDistance(slice, front_points_2d(1, :), front_points_2d(3, :));

% figure;
% plot(slice(:, 1), slice(:, 2));
% hold on;
% plot(front_points_2d(:, 1), front_points_2d(:, 2), 'kx');
% hold off;

end
