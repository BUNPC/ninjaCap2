function [v_new, e_new, vOut_new, eOut_new] = reclassifyInternalBoundaryEdges(v, e, vOut, eOut, polygon_coords)
%reclassifyInternalBoundaryEdges Moves internal edges from boundary data to main mesh data.
%
%   INPUTS:
%   v, e           - The main mesh vertices and edges.
%   vOut, eOut     - The boundary data vertices and edges.
%   polygon_coords - Px2 matrix of [x, y] coordinates defining the polygon.
%
%   OUTPUTS:
%   v_new, e_new     - The updated main mesh with the reclassified edges/vertices added.
%   vOut_new, eOut_new - The updated boundary data with the internal edges/vertices removed.

PROXIMITY_TOLERANCE = 1e-4;

% Initialize new variables with the originals
v_new = v;
e_new = e;

%% Step 1: Identify which eOut edges are fully internal
if isempty(eOut) || isempty(vOut)
    % No boundary edges to process, return the inputs as is
    vOut_new = vOut;
    eOut_new = eOut;
    return;
end

is_inside = inpolygon(vOut(:,1), vOut(:,2), polygon_coords(:,1), polygon_coords(:,2));
v1_is_inside = is_inside(eOut(:,1));
v2_is_inside = is_inside(eOut(:,2));
is_edge_internal = v1_is_inside & v2_is_inside;

%% Step 2: Partition boundary data into "to move" and "to keep"
eOut_to_move = eOut(is_edge_internal, :);
eOut_to_keep = eOut(~is_edge_internal, :);

% Get the original indices and coordinates of vertices that need to be moved
v_indices_to_move = unique(eOut_to_move(:));
v_coords_to_move = vOut(v_indices_to_move, :);

%% Step 3: Add the "moved" vertices and edges to the main mesh
map_vout_to_vnew = zeros(size(vOut, 1), 1);

% Merge the vertices, checking for duplicates
for i = 1:length(v_indices_to_move)
    old_vout_idx = v_indices_to_move(i);
    current_coord = vOut(old_vout_idx, :);
    
    % Check if this vertex coordinate already exists in the main vertex list
    if ~isempty(v_new)
        dists = vecnorm(v_new - current_coord, 2, 2);
        [min_dist, existing_idx] = min(dists);
    else
        min_dist = inf;
        existing_idx = 0;
    end
    
    if min_dist < PROXIMITY_TOLERANCE
        % It already exists, map to the existing index
        map_vout_to_vnew(old_vout_idx) = existing_idx;
    else
        % It's a new vertex, append it and map to the new index
        v_new(end+1, :) = current_coord;
        map_vout_to_vnew(old_vout_idx) = size(v_new, 1);
    end
end

% Re-index the edges that are being moved and add them to the main edge list
if ~isempty(eOut_to_move)
    e_to_add = map_vout_to_vnew(eOut_to_move);
    e_new = [e_new; e_to_add];
end

%% Step 4: Clean up the remaining boundary data
% The new boundary edge list is just the edges we decided to keep
eOut_new_raw = eOut_to_keep;

% Now, clean up any orphan nodes in vOut left by the removal of edges
if isempty(eOut_new_raw)
    vOut_new = [];
    eOut_new = [];
else
    used_vout_indices = unique(eOut_new_raw(:));
    vOut_new = vOut(used_vout_indices, :);
    
    map_old_vout_to_new_vout = zeros(size(vOut, 1), 1);
    map_old_vout_to_new_vout(used_vout_indices) = 1:length(used_vout_indices);
    
    eOut_new = map_old_vout_to_new_vout(eOut_new_raw);
end

end