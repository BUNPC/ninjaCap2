function o = createHybridCapDesign_v2(sideLeftOutline, sideLeftIdx, sideRightOutline, sideRightIdx, topOutline, topIdx)
% createHybridCapDesign_v2 - Fills cap outlines with a hybrid geometry.
%
% This version creates exclusion zones around the optode geometry BEFORE
% filling the remaining space with a hexagonal lattice.

%% --- Initial Setup & Parameters ---
o(1).v = [50 50]; o(1).xoutline = topOutline(:,1) - min(topOutline(:,1)); o(1).youtline = topOutline(:,2) - min(topOutline(:,2));
o(2).v = [100 100]; o(2).xoutline = sideLeftOutline(:,1) - min(sideLeftOutline(:,1)); o(2).youtline = sideLeftOutline(:,2) - min(sideLeftOutline(:,2));
o(3).v = [100 100]; o(3).xoutline = sideRightOutline(:,1) - min(sideRightOutline(:,1)); o(3).youtline = sideRightOutline(:,2) - min(sideRightOutline(:,2));

hEdge = 10; % Hexagon edge length
nOutlines = length(o);
clearance = hEdge; % Define keep-out distance around optode structures

%% --- Generate Priority Optode Geometry ---
% Call the function to get all optode-related points and connections
[~, sources_left, ~, sources_right, ...
 det_centroids_left, src_centroids_left, det_centroids_right, src_centroids_right, ...
 outer_sources_left_side, outer_sources_right_side, outer_detectors_left_side, outer_detectors_right_side] = ...
 get_UHD_optode_pos_modified(sideLeftOutline(sideLeftIdx(1):sideLeftIdx(2),:), ...
                             sideRightOutline(sideRightIdx(1):sideRightIdx(2),:));

% Generate connection lines
conn_det_left = find_connections(det_centroids_left, sources_left, 3, 15);
conn_src_left = find_connections(src_centroids_left, sources_left, 3, 15);
all_connections_left_raw = [conn_det_left; conn_src_left];

conn_det_right = find_connections(det_centroids_right, sources_right, 3, 15);
conn_src_right = find_connections(src_centroids_right, sources_right, 3, 15);
all_connections_right_raw = [conn_det_right; conn_src_right];

% Define all points to be included in the exclusion zone
exclusion_pts_left = [outer_sources_left_side; outer_detectors_left_side; det_centroids_left; src_centroids_left];
exclusion_pts_right = [outer_sources_right_side; outer_detectors_right_side; det_centroids_right; src_centroids_right];

% Convert connections to vertex/edge format and shift to local coordinates
offset_left = min(sideLeftOutline(:,1:2));
[v_conn_left, e_conn_left] = convert_plot_lines_to_graph(all_connections_left_raw);
v_conn_left = v_conn_left - offset_left;
exclusion_pts_left = exclusion_pts_left - offset_left;
all_connections_left_shifted = all_connections_left_raw;
if ~isempty(all_connections_left_shifted)
    all_connections_left_shifted(~isnan(all_connections_left_shifted(:,1)),:) = all_connections_left_shifted(~isnan(all_connections_left_shifted(:,1)),:) - offset_left;
end

offset_right = min(sideRightOutline(:,1:2));
[v_conn_right, e_conn_right] = convert_plot_lines_to_graph(all_connections_right_raw);
v_conn_right = v_conn_right - offset_right;
exclusion_pts_right = exclusion_pts_right - offset_right;
all_connections_right_shifted = all_connections_right_raw;
if ~isempty(all_connections_right_shifted)
    all_connections_right_shifted(~isnan(all_connections_right_shifted(:,1)),:) = all_connections_right_shifted(~isnan(all_connections_right_shifted(:,1)),:) - offset_right;
end


%% --- Set up Seams (Unchanged) ---
% This entire section is identical to the previous version
iSeam = 1;
xx = o(iSeam).xoutline; yy = o(iSeam).youtline;
o(iSeam).seamLength = zeros(length(topOutline),1); o(iSeam).seamIs = zeros(length(topOutline),1); o(iSeam).seamIs(topIdx(1)) = 1;
for ii = topIdx(1)+1:topIdx(2)
    o(iSeam).seamLength(ii) = o(iSeam).seamLength(ii-1) + norm( [xx(ii)-xx(ii-1) yy(ii)-yy(ii-1)] ); o(iSeam).seamIs(ii) = 1;
end
o(iSeam).seamIs(topIdx(3)) = 1;
for ii = topIdx(3)+1:topIdx(4)
    o(iSeam).seamLength(ii) = o(iSeam).seamLength(ii-1) + norm( [xx(ii)-xx(ii-1) yy(ii)-yy(ii-1)] ); o(iSeam).seamIs(ii) = 1;
end
iSeam = 2;
o(iSeam).seamStartIdx = sideLeftIdx(1); o(iSeam).seamEndIdx = sideLeftIdx(2); o(iSeam).seamLength = zeros(length(sideLeftOutline),1); o(iSeam).seamIs = zeros(length(sideLeftOutline),1);
xx = o(iSeam).xoutline; yy = o(iSeam).youtline; o(iSeam).seamIs(o(iSeam).seamStartIdx) = 1;
for ii = o(iSeam).seamStartIdx+1 : o(iSeam).seamEndIdx
    o(iSeam).seamLength(ii) = o(iSeam).seamLength(ii-1) + norm( [xx(ii)-xx(ii-1) yy(ii)-yy(ii-1)] ); o(iSeam).seamIs(ii) = 1;
end
iSeam = 3;
o(iSeam).seamStartIdx = sideRightIdx(1); o(iSeam).seamEndIdx = sideRightIdx(2); o(iSeam).seamLength = zeros(length(sideRightOutline),1); o(iSeam).seamIs = zeros(length(sideRightOutline),1);
xx = o(iSeam).xoutline; yy = o(iSeam).youtline; o(iSeam).seamIs(o(iSeam).seamStartIdx) = 1;
for ii = o(iSeam).seamStartIdx+1 : o(iSeam).seamEndIdx
    o(iSeam).seamLength(ii) = o(iSeam).seamLength(ii-1) + norm( [xx(ii)-xx(ii-1) yy(ii)-yy(ii-1)] ); o(iSeam).seamIs(ii) = 1;
end


%% --- NEW: Create Mask, Apply Exclusions, and Fill with Hexagons ---
for iO = 1:nOutlines
    % Create the base mask from the outline
    [xgrid,ygrid] = meshgrid(0:max(o(iO).xoutline)*1.05, 0:max(o(iO).youtline)*1.05);
    base_mask = inpolygon(xgrid,ygrid,o(iO).xoutline,o(iO).youtline);
    
    final_mask = base_mask;
    
    % For side panels, modify the mask to exclude optode zones
    if iO == 2 % Left Panel
        final_mask = create_exclusion_mask(base_mask, xgrid, ygrid, exclusion_pts_left, all_connections_left_shifted, clearance);
        v_conn = v_conn_left;
        e_conn = e_conn_left;
    elseif iO == 3 % Right Panel
        final_mask = create_exclusion_mask(base_mask, xgrid, ygrid, exclusion_pts_right, all_connections_right_shifted, clearance);
        v_conn = v_conn_right;
        e_conn = e_conn_right;
    end
    
    % Fill the final, possibly modified, mask with hexagons
    [v_hex, e_hex, vOut_hex, eOut_hex] = fillCapWithHexagons_func(o(iO).v, hEdge, final_mask);
    
    if iO == 1 % Top Panel
        o(iO).v = v_hex;
        o(iO).e = e_hex;
    else % Side Panels: Combine the hex fill and the connection geometry
        num_hex_verts = size(v_hex, 1);
        if ~isempty(e_conn)
            e_conn_remapped = e_conn + num_hex_verts;
            o(iO).v = [v_hex; v_conn];
            o(iO).e = [e_hex; e_conn_remapped];
        else
            o(iO).v = v_hex;
            o(iO).e = e_hex;
        end
    end
    
    o(iO).vOut = vOut_hex;
    o(iO).eOut = eOut_hex;
    o(iO).Imask = final_mask; % Store the final mask for visualization if needed
end

%% --- Stitching and Cleanup (Unchanged) ---
% The rest of the script (pulling vertices to outline, crossing seams,
% removing duplicates, etc.) is identical to the original and will now
% operate on the new hybrid geometry.

% (The rest of your original code for stitching and cleanup follows here)
% ...

end % End of main function

% -------------------------------------------------------------------------
% LOCAL HELPER FUNCTIONS
% -------------------------------------------------------------------------

function modified_mask = create_exclusion_mask(base_mask, xgrid, ygrid, ex_pts, ex_lines, clearance)
% Removes areas from a mask that are too close to specified points and lines.

    if isempty(ex_pts) && isempty(ex_lines)
        modified_mask = base_mask;
        return;
    end

    % Get coordinates of grid points inside the initial mask
    [rows, cols] = find(base_mask);
    if isempty(rows)
        modified_mask = base_mask;
        return;
    end
    
    % Efficiently get the unique x and y coordinates of the grid points
    grid_pts_x = xgrid(1, cols)';
    grid_pts_y = ygrid(rows, 1);
    % Create the list of points to check
    grid_pts = [grid_pts_x, grid_pts_y(1:length(grid_pts_x))]; % Minor fix for non-square grids
    
    to_remove_mask = false(size(grid_pts, 1), 1);

    % 1. Find grid points to remove based on proximity to exclusion points
    if ~isempty(ex_pts)
        min_dist_to_points = min(pdist2(grid_pts, ex_pts), [], 2);
        to_remove_mask = to_remove_mask | (min_dist_to_points < clearance);
    end
    
    % 2. Find grid points to remove based on proximity to exclusion lines
    if ~isempty(ex_lines)
        conn_segments = [];
        for i = 1:3:size(ex_lines, 1)
            segment = ex_lines(i:i+1, :);
            if ~any(isnan(segment), 'all')
                conn_segments(end+1, :, :) = segment;
            end
        end
        
        if ~isempty(conn_segments)
            min_dist_to_lines = inf(size(grid_pts, 1), 1);
            for i = 1:size(conn_segments, 1)
                v1 = squeeze(conn_segments(i, 1, :))';
                v2 = squeeze(conn_segments(i, 2, :))';
                dists_to_segment = point_to_segment_distance_vectorized(grid_pts, v1, v2);
                min_dist_to_lines = min(min_dist_to_lines, dists_to_segment);
            end
            to_remove_mask = to_remove_mask | (min_dist_to_lines < clearance);
        end
    end
    
    % Update the final mask
    modified_mask = base_mask;
    linear_indices_to_remove = sub2ind(size(base_mask), rows(to_remove_mask), cols(to_remove_mask));
    modified_mask(linear_indices_to_remove) = 0;
end

function dists = point_to_segment_distance_vectorized(pts, v1, v2)
% Vectorized calculation of shortest distance from points to a line segment.
    v1v2 = v2 - v1;
    ptv1 = pts - v1;
    L2 = sum(v1v2.^2);
    if L2 == 0, dists = sqrt(sum(ptv1.^2, 2)); return; end
    t = (ptv1 * v1v2') / L2;
    t = max(0, min(1, t));
    projections = v1 + t .* v1v2;
    dists = sqrt(sum((pts - projections).^2, 2));
end



% --- Also include your other original helper functions ---
% find_connections, get_det_centroids, point_to_segment_distance, etc.
% convert_plot_lines_to_graph

function [v, e] = convert_plot_lines_to_graph(lines)
% Converts a list of line segments for plotting into a vertex and edge list
    if isempty(lines)
        v = []; e = []; return;
    end
    valid_lines = lines(~isnan(lines(:,1)), :);
    if isempty(valid_lines)
        v = []; e = []; return;
    end
    
    starts = valid_lines(1:2:end, :);
    ends = valid_lines(2:2:end, :);
    all_points = [starts; ends];
    
    [v, ~, ic] = unique(all_points, 'rows', 'stable');
    
    num_starts = size(starts, 1);
    e = [ic(1:num_starts), ic(num_starts+1:end)];
end