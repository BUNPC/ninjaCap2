function [v_trimmed, e_trimmed] = trimMeshToPolygon(v, e, polygon_coords)
%trimMeshToPolygon Removes edges and vertices that are fully outside a polygon.
%
%   This function takes a mesh (defined by vertices 'v' and edges 'e') and
%   a polygon. It removes any edge from 'e' where BOTH of its vertices
%   fall outside the polygon boundary. After removing the external edges,
%   it also removes any vertices from 'v' that are no longer referenced
%   by any edge.
%
%   INPUTS:
%   v              - Nx2 matrix of [x, y] vertex coordinates.
%   e              - Mx2 matrix of edges, with indices corresponding to 'v'.
%   polygon_coords - Px2 matrix of [x, y] coordinates defining the polygon.
%
%   OUTPUTS:
%   v_trimmed      - The final list of trimmed vertex coordinates.
%   e_trimmed      - The final, re-indexed list of trimmed edges.

% If there are no edges to start with, return empty lists.
if isempty(e)
    v_trimmed = v; % Or [] depending on desired behavior
    e_trimmed = [];
    return;
end

%% Step 1: Identify Vertices Outside the Polygon
% Use MATLAB's inpolygon function to determine which vertices are inside.
% 'is_inside' will be a logical vector (true for inside, false for outside).
is_inside = inpolygon(v(:,1), v(:,2), polygon_coords(:,1), polygon_coords(:,2));
is_outside = ~is_inside;
outside_indices = find(is_outside);

%% Step 2: Identify Edges to Remove
% An edge should be removed if its start vertex AND end vertex are outside.
v1_is_outside = is_outside(e(:,1));
v2_is_outside = is_outside(e(:,2));

is_edge_to_remove = v1_is_outside & v2_is_outside;

% Keep all edges that are NOT marked for removal.
e_kept = e(~is_edge_to_remove, :);

%% Step 3: Remove Orphan Vertices and Re-index
% After removing edges, some vertices may no longer be part of the mesh.
% We must remove these "orphan" vertices and update the edge indices.

% Find all unique vertex indices that are still used in the kept edges.
used_indices = unique(e_kept(:));

% Create the final, trimmed vertex list.
v_trimmed = v(used_indices, :);

% To re-index the edges, create a map from the old vertex indices to the new ones.
map_old_to_new = zeros(size(v, 1), 1);
map_old_to_new(used_indices) = 1:length(used_indices);

% Apply the map to the kept edge list to get the final, trimmed edge list.
e_trimmed = map_old_to_new(e_kept);

% Final cleanup for consistent ordering.
e_trimmed = unique(sort(e_trimmed, 2), 'rows');

end