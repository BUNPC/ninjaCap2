function new_plane = planeMoveToIntersect(plane, pt)
%PLANEMOVETOINTERSECT Summary of this function goes here
%   Detailed explanation goes here

% plane: [x0 y0 z0 dx1 dy1 dz1 dx2 dy2 dz2]
% based on: matGeom implementation
plane_origin = plane(1:3);
plane_normal = cross(plane(4:6), plane(7:9));
plane_normal = plane_normal ./ norm(plane_normal); % normal

% reference:
% https://stackoverflow.com/questions/9605556/how-to-project-a-point-onto-a-plane-in-3d#9605695

% distance of vertices from plane along normal
dist = dot(plane_normal, pt - plane_origin, 2);

% shift origin
new_plane = plane + [(dist * plane_normal) 0 0 0 0 0 0];

end

