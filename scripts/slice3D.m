function edges = slice3D(v, f, plane)
%SLICE3D Summary of this function goes here
%   Detailed explanation goes here

% plane: [x0 y0 z0 dx1 dy1 dz1 dx2 dy2 dz2]
% based on: matGeom implementation
plane_origin = plane(1:3);
plane_basis1 = plane(4:6) ./ norm(plane(4:6));
plane_basis2 = plane(7:9) ./ norm(plane(7:9));
plane_normal = cross(plane_basis1, plane_basis2);
plane_normal = plane_normal ./ norm(plane_normal);

% algorithm based on: http://geomalgorithms.com/a06-_intersect-2.html#Triangle-Plane

% position of vertices relative to the plane
s = dot(repmat(plane_normal, size(v, 1), 1), v - plane_origin, 2);

% get faces of interest
f_s = s(f);
faces_spanning = and(any(f_s <= 0, 2), any(f_s >= 0, 2));

% edges [x1, y1, z1; z2, y2, z2]
edges_3d = [];

% for each spanning face
for i = find(faces_spanning)'
    % directions are relative, just helpful to have names
    left = f_s(i, :) < 0;
    right = f_s(i, :) > 0;
    
    num_left = sum(left);
    num_right = sum(right);
    num_on = 3 - num_left - num_right;
    
    if num_on == 0
        % most common scenario
        
        % find vertices of two lines that intersect plane
        idx = find(left([1 2 3]) ~= left([2 3 1]));
        assert(length(idx) == 2);
        idx_line1_point1 = 1 + mod(idx(1) - 1, 3);
        idx_line1_point2 = 1 + mod(idx(1), 3);
        idx_line2_point1 = 1 + mod(idx(2) - 1, 3);
        idx_line2_point2 = 1 + mod(idx(2), 3);
        
        line1_origin = v(f(i, idx_line1_point1), :);
        line1_slope = v(f(i, idx_line1_point2), :) - line1_origin;
        line1_r = dot(plane_normal, plane_origin - line1_origin) ./ dot(plane_normal, line1_slope);
        point1 = line1_origin + line1_r * line1_slope;
        
        line2_origin = v(f(i, idx_line2_point1), :);
        line2_slope = v(f(i, idx_line2_point2), :) - line2_origin;
        line2_r = dot(plane_normal, plane_origin - line2_origin) ./ dot(plane_normal, line2_slope);
        point2 = line2_origin + line2_r * line2_slope;
        
        edges_3d = [edges_3d; point1; point2];
    elseif num_on == 1
        % only add if other two points spans, othwerwise intersection is
        % just a single point
        if num_left == 1 && num_right == 1
            line1_origin = v(f(i, left), :);
            line1_slope = v(f(i, right), :) - line1_origin;
            line1_r = dot(plane_normal, plane_origin - line1_origin) ./ dot(plane_normal, line1_slope);
            point1 = line1_origin + line1_r * line1_slope;
            
            point2 = v(f(i, ~left & ~right), :);
            
            edges_3d = [edges_3d; point1; point2];
        end
    elseif num_on == 2
        % edge between two points on the plane
        idx = find(~left & ~right);
        edges_3d = [edges_3d; v(f(i, idx(1))); v(f(i, idx(2)))];
    elseif num_on == 3
        % add all three edges
        edges_3d = [edges_3d; v(f(i, 1), :); v(f(i, 2), :); v(f(i, 2), :); v(f(i, 3), :); v(f(i, 3), :); v(f(i, 1), :)];
    end
end

% project to 2d
% based on: https://stackoverflow.com/questions/23472048/projecting-3d-points-to-2d-plane#23474396

edges_dim1 = dot(repmat(plane_basis1, size(edges_3d, 1), 1), edges_3d - plane_origin, 2);
edges_dim2 = dot(repmat(plane_basis2, size(edges_3d, 1), 1), edges_3d - plane_origin, 2);

edges = nan(size(edges_3d, 1) * 3 / 2, 2);
edges(1:3:end, 1) = edges_dim1(1:2:end);
edges(2:3:end, 1) = edges_dim1(2:2:end);
edges(1:3:end, 2) = edges_dim2(1:2:end);
edges(2:3:end, 2) = edges_dim2(2:2:end);

end

