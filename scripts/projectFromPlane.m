function pts = projectFromPlane(plane, pts_projected_2d)
%PROJECTFROMPLANE Summary of this function goes here
%   Detailed explanation goes here

% plane: [x0 y0 z0 dx1 dy1 dz1 dx2 dy2 dz2]
% based on: matGeom implementation
plane_origin = plane(1:3);
plane_basis1 = plane(4:6) ./ norm(plane(4:6));
plane_basis2 = plane(7:9) ./ norm(plane(7:9));
%plane_normal = cross(plane_basis1, plane_basis2);
%plane_normal = plane_normal ./ norm(plane_normal);

% reference:
% https://stackoverflow.com/questions/9605556/how-to-project-a-point-onto-a-plane-in-3d#9605695

% reverse 2d projection
pts = bsxfun(@plus, plane_origin, pts_projected_2d(:, 1) * plane_basis1 + pts_projected_2d(:, 2) * plane_basis2);

end
