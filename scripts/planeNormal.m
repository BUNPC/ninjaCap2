function plane_normal = planeNormal(plane)
%PLANENORMAL Summary of this function goes here
%   Detailed explanation goes here

% plane: [x0 y0 z0 dx1 dy1 dz1 dx2 dy2 dz2]
% based on: matGeom implementation
%plane_origin = plane(1:3);
plane_normal = cross(plane(4:6), plane(7:9));
plane_normal = plane_normal ./ norm(plane_normal); % normal

end

