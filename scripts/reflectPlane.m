function pts_reflected = reflectPlane(plane, pts)
%PLANEREFLECT Summary of this function goes here
%   Detailed explanation goes here

% plane: [x0 y0 z0 dx1 dy1 dz1 dx2 dy2 dz2]
% based on: matGeom implementation
plane_origin = plane(1:3);
plane_normal = cross(plane(4:6), plane(7:9));
plane_normal = plane_normal ./ norm(plane_normal); % normal

% reference:
% https://stackoverflow.com/questions/9605556/how-to-project-a-point-onto-a-plane-in-3d#9605695

% distance of vertices from plane along normal
dist = dot(repmat(plane_normal, size(pts, 1), 1), pts - plane_origin, 2);

% reflect over plane (1 * dist * plane_normal is distance to plane)
pts_reflected = pts - (2 .* dist .* plane_normal);

end
