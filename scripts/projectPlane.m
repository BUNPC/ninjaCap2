function [pts_projected, pts_projected_2d] = projectPlane(plane, pts)
%PROJECTPLANE Summary of this function goes here
%   Detailed explanation goes here

% plane: [x0 y0 z0 dx1 dy1 dz1 dx2 dy2 dz2]
% based on: matGeom implementation
plane_origin = plane(1:3);
plane_basis1 = plane(4:6) ./ norm(plane(4:6));
plane_basis2 = plane(7:9) ./ norm(plane(7:9));
plane_normal = cross(plane_basis1, plane_basis2);
plane_normal = plane_normal ./ norm(plane_normal);

% reference:
% https://stackoverflow.com/questions/9605556/how-to-project-a-point-onto-a-plane-in-3d#9605695

% distance of vertices from plane along normal
dist = dot(repmat(plane_normal, size(pts, 1), 1), pts - plane_origin, 2);

% project onto plane
pts_projected = pts - dist .* plane_normal;

% 2d projection
if nargout > 1
    pts_projected_2d = zeros(size(pts_projected, 1), 2);
    pts_projected_2d(:, 1) = dot(repmat(plane_basis1, size(pts_projected, 1), 1), pts_projected - plane_origin, 2);
    pts_projected_2d(:, 2) = dot(repmat(plane_basis2, size(pts_projected, 1), 1), pts_projected - plane_origin, 2);
end

end
