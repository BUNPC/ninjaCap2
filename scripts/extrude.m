function [faces, vertices] = extrude(poly, height)
%EXTRUDE Summary of this function goes here
%   Detailed explanation goes here

% figure out triangulation
tri = triangulation(poly);

% get points and faces
vertices = tri.Points;
faces = tri.ConnectivityList;

% convert to 3D
vertices = [vertices zeros(size(vertices, 1), 1)];

%% should point outside (negative)

% calculate normals
normals = cross(vertices(faces(:, 2), :) - vertices(faces(:, 1), :), vertices(faces(:, 3), :) - vertices(faces(:, 1), :));

% flip faces
to_flip = normals(:, 3) < 0;
faces(to_flip, :) = fliplr(faces(to_flip, :));

%% build walls

if height == 0
    return;
end

% look up for first set of vertices (z = 0)
[found, v_idx] = ismember(poly.Vertices, tri.Points, 'rows');
found(end + 1) = false; % extra one for convenience

% second set of vertices (z = height)
num_vertices = size(vertices, 1);
vertices = [vertices; bsxfun(@plus, vertices, [0 0 height])];

% second set of faces connecting extruded vertices
faces = [faces; fliplr(faces + num_vertices)];

last_start = 1;
for i = 1:length(found)
    if found(i)
        f1 = v_idx(i);
        if found(i + 1)
            f2 = v_idx(i + 1);
        else
            f2 = v_idx(last_start);
        end
        
        % triangles for wall
        t = [f1 f1 + num_vertices f2; ...
            f1 + num_vertices f2 + num_vertices f2];
        
        % normal
        n = cross(vertices(t(:, 2), :) - vertices(t(:, 1), :), vertices(t(:, 3), :) - vertices(t(:, 1), :));
        
        % check normal (generate midpoints and move along normal)
        p = (vertices(f1, 1:2) + vertices(f2, 1:2)) ./ 2;
        p = p + n(1, 1:2) * eps;
        
        % flip if pointing inward
        if ~isinterior(poly, p)
            t = fliplr(t);
        end
        
        % add faces
        faces = [faces; t];
    else
        last_start = i + 1;
    end
end

end
