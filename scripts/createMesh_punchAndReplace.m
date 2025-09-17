function [v_final, e_final, vOut_final, eOut_final] = createMesh_punchAndReplace(Imask, hEdge, v_conn, e_conn, punch_radius)
%createMesh_punchAndReplace that trims the predefined structure at the boundary.

PROXIMITY_TOLERANCE = 1e-4;

%% Task 1: Full Hexagonal Fill (includes vOut/eOut calculation)
disp('Task 1: Filling mask and calculating initial boundary data...');
v_full = []; e_full = []; vOut_full = []; eOut_full = [];
[valid_y, valid_x] = find(Imask);
if isempty(valid_y), v_final=v_conn; e_final=e_conn; vOut_final=[]; eOut_final=[]; return; end
v_seed = [valid_x(round(end/2)), valid_y(round(end/2))];
v_full(1,:) = v_seed; vlst = [1]; voe(1,1) = 0;
[ny, nx] = size(Imask);
dv = [hEdge 0; hEdge*cos(2*pi/3) hEdge*sin(2*pi/3); hEdge*cos(2*pi/3) -hEdge*sin(2*pi/3)];
while ~isempty(vlst)
    current_v_idx = vlst(end); vlst(end) = [];
    if voe(current_v_idx) == 0, vtmp = v_full(current_v_idx, :) + dv; else, vtmp = v_full(current_v_idx, :) - dv; end
    oe = mod(voe(current_v_idx) + 1, 2);
    for ii = 1:size(vtmp, 1)
        p_new = vtmp(ii, :);
        dists = vecnorm(v_full - p_new, 2, 2);
        [min_dist, closest_idx] = min(dists);
        is_inside = p_new(1)>=1 && p_new(1)<=nx && p_new(2)>=1 && p_new(2)<=ny && Imask(round(p_new(2)), round(p_new(1)))==1;
        if min_dist < PROXIMITY_TOLERANCE
            if isempty(e_full) || ~ismember(sort([current_v_idx, closest_idx]), sort(e_full, 2), 'rows'), e_full(end+1, :) = [current_v_idx, closest_idx]; end
        elseif is_inside
            v_full(end+1, :) = p_new; new_v_idx = size(v_full, 1);
            e_full(end+1, :) = [current_v_idx, new_v_idx];
            voe(new_v_idx) = oe; vlst(end+1) = new_v_idx;
        else
            p1 = v_full(current_v_idx, :); p2 = p_new;
            if ~isempty(vOut_full), ii1 = find(vecnorm(vOut_full-p1,2,2) < PROXIMITY_TOLERANCE, 1); else, ii1 = []; end
            if isempty(ii1), vOut_full(end+1,:) = p1; ii1 = size(vOut_full,1); end
            if ~isempty(vOut_full), ii2 = find(vecnorm(vOut_full-p2,2,2) < PROXIMITY_TOLERANCE, 1); else, ii2 = []; end
            if isempty(ii2), vOut_full(end+1,:) = p2; ii2 = size(vOut_full,1); end
            if isempty(eOut_full) || ~ismember(sort([ii1 ii2]), sort(eOut_full,2), 'rows'), eOut_full(end+1,:) = [ii1 ii2]; end
        end
    end
end

%% Task 2: Punch the Hole in the Grid
disp('Task 2: Punching hole for new structure...');
punch_mask = createExclusionMask(v_conn, e_conn, size(Imask), punch_radius);
v_full_round = round(v_full);
v_full_round(:,1) = max(1, min(v_full_round(:,1), size(Imask,2)));
v_full_round(:,2) = max(1, min(v_full_round(:,2), size(Imask,1)));
linear_indices = sub2ind(size(Imask), v_full_round(:,2), v_full_round(:,1));
is_vertex_to_remove = punch_mask(linear_indices);
v_indices_to_remove = find(is_vertex_to_remove);
is_edge_to_remove = ismember(e_full(:,1), v_indices_to_remove) | ismember(e_full(:,2), v_indices_to_remove);
all_edge_verts = e_full(is_edge_to_remove, :);
boundary_indices_old = setdiff(unique(all_edge_verts(:)), v_indices_to_remove);
v_indices_to_keep = find(~is_vertex_to_remove);
v_punched = v_full(v_indices_to_keep, :);
e_to_remap = e_full(~is_edge_to_remove, :);
map_old_to_new = zeros(size(v_full, 1), 1);
map_old_to_new(v_indices_to_keep) = 1:length(v_indices_to_keep);
e_punched = map_old_to_new(e_to_remap);
boundary_indices_new = map_old_to_new(boundary_indices_old);

%% Task 3: Insert and Stitch the Meshes
disp('Task 3: Stitching meshes with zipper method...');
num_punched_verts = size(v_punched, 1);
v_temp = [v_punched; v_conn];
e_temp = e_punched; % NOTE: e_conn is NOT added yet
stitching_edges = [];
available_hex_indices = boundary_indices_new;
available_vconn_local_indices = (1:size(v_conn, 1))';
while ~isempty(available_hex_indices) && ~isempty(available_vconn_local_indices)
    available_hex_coords = v_punched(available_hex_indices, :);
    available_vconn_coords = v_conn(available_vconn_local_indices, :);
    distances = pdist2(available_hex_coords, available_vconn_coords);
    [~, min_linear_idx] = min(distances(:));
    [hex_list_idx, vconn_list_idx] = ind2sub(size(distances), min_linear_idx);
    chosen_hex_idx = available_hex_indices(hex_list_idx);
    chosen_vconn_local_idx = available_vconn_local_indices(vconn_list_idx);
    chosen_vconn_final_idx = num_punched_verts + chosen_vconn_local_idx;
    stitching_edges(end+1, :) = [chosen_hex_idx, chosen_vconn_final_idx];
    available_hex_indices(hex_list_idx) = [];
    available_vconn_local_indices(vconn_list_idx) = [];
end
e_temp = [e_temp; stitching_edges];

%% Task 4: Process Predefined Structure and Finalize Boundary Data
disp('Task 4: Trimming predefined structure and finalizing boundary data...');

Imask_double = double(Imask);
e_conn_inside = [];
vOut_conn_temp = []; eOut_conn_temp = [];

% Part A: Partition e_conn into "inside" and "outside" sets
for k = 1:size(e_conn, 1)
    p1_idx = e_conn(k, 1); p2_idx = e_conn(k, 2);
    p1 = v_conn(p1_idx, :); p2 = v_conn(p2_idx, :);
    line_pixels = improfile(Imask_double, [p1(1), p2(1)], [p1(2), p2(2)]);
    
    if any(line_pixels == 0) || any(isnan(line_pixels))
        % This is an OUTSIDE edge
        if ~isempty(vOut_conn_temp), ii1=find(vecnorm(vOut_conn_temp-p1,2,2)<PROXIMITY_TOLERANCE,1); else, ii1=[]; end
        if isempty(ii1), vOut_conn_temp(end+1,:)=p1; ii1=size(vOut_conn_temp,1); end
        if ~isempty(vOut_conn_temp), ii2=find(vecnorm(vOut_conn_temp-p2,2,2)<PROXIMITY_TOLERANCE,1); else, ii2=[]; end
        if isempty(ii2), vOut_conn_temp(end+1,:)=p2; ii2=size(vOut_conn_temp,1); end
        eOut_conn_temp(end+1,:) = [ii1 ii2];
    else
        % This is an INSIDE edge
        e_conn_inside(end+1, :) = e_conn(k, :);
    end
end
% Add the "inside" edges to the main temporary edge list
e_temp = [e_temp; e_conn_inside + num_punched_verts];

% Part B: Filter the original hex boundary data
removed_hex_coords = v_full(v_indices_to_remove, :);
is_vout_to_remove = false(size(vOut_full, 1), 1);
if ~isempty(removed_hex_coords) && ~isempty(vOut_full), dists = pdist2(vOut_full, removed_hex_coords); is_vout_to_remove = min(dists,[],2) < PROXIMITY_TOLERANCE; end
vOut_indices_to_remove = find(is_vout_to_remove);
is_eout_to_remove = ismember(eOut_full(:,1), vOut_indices_to_remove) | ismember(eOut_full(:,2), vOut_indices_to_remove);
eOut_hex_kept_remapped = eOut_full(~is_eout_to_remove, :);
vOut_indices_to_keep = find(~is_vout_to_remove);
vOut_hex_kept = vOut_full(vOut_indices_to_keep, :);
map_old_vout_to_new = zeros(size(vOut_full, 1), 1);
map_old_vout_to_new(vOut_indices_to_keep) = 1:length(vOut_indices_to_keep);
eOut_hex_kept = map_old_vout_to_new(eOut_hex_kept_remapped);

% Part C: Merge the two boundary data sets
vOut_final = vOut_hex_kept;
eOut_final = eOut_hex_kept;
map_conn_to_final = zeros(size(vOut_conn_temp, 1), 1);
for i = 1:size(vOut_conn_temp, 1)
    pt = vOut_conn_temp(i, :);
    final_idx = find(vecnorm(vOut_final - pt, 2, 2) < PROXIMITY_TOLERANCE, 1);
    if isempty(final_idx), vOut_final(end+1,:) = pt; final_idx = size(vOut_final, 1); end
    map_conn_to_final(i) = final_idx;
end
if ~isempty(eOut_conn_temp)
    eOut_conn_final = map_conn_to_final(eOut_conn_temp);
    eOut_final = [eOut_final; eOut_conn_final];
end

%% Task 5: Final Cleanup of Orphan Nodes
disp('Task 5: Cleaning up unreferenced nodes...');
used_indices = unique(e_temp(:));
v_final = v_temp(used_indices, :);
map_temp_to_final = zeros(size(v_temp, 1), 1);
map_temp_to_final(used_indices) = 1:length(used_indices);
e_final = map_temp_to_final(e_temp);

% Final cleanup
if ~isempty(e_final), e_final = unique(sort(e_final, 2), 'rows'); end
if ~isempty(eOut_final), eOut_final = unique(sort(eOut_final, 2), 'rows'); end
disp('Process complete.');
end