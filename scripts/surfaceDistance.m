function distance = surfaceDistance(vertices, faces, p1, p2)
%SURFACEDISTANCE Summary of this function goes here
%   Detailed explanation goes here

% get mid point
p_mid = closestVertex(vertices, (p1 + p2) ./ 2);
p1 = closestVertex(vertices, p1);
p2 = closestVertex(vertices, p2);

% plane
plane = makePlane([p1; p_mid; p2]);

% slice
slice = slice3D(vertices, faces, plane);

% projected
[~, proj_point] = projectPlane(plane, [p1; p2]);

% % debug
% figure;
% plot(slice(:, 1), slice(:, 2));
% hold on;
% scatter(proj_point(:, 1), proj_point(:, 2), 'filled');
% hold off;

% calculate distance
distance = sliceDistance(slice, proj_point(1, :), proj_point(2, :));

end

