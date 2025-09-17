function o = createHybridCapDesign(sideLeftOutline, sideLeftIdx, sideRightOutline, sideRightIdx, topOutline, topIdx)
% This version includes debugging plots to help isolate visual errors.

%% --- Configuration ---
hEdge = 10; % Hexagon edge length
clearance = hEdge / 4; % Keep-out distance around optode structures
connection_threshold = hEdge * 1.5; % Max distance to add a stitching edge

%% --- Initial Struct Setup ---
o(1).v = [50 50]; o(1).xoutline = topOutline(:,1) - min(topOutline(:,1)); o(1).youtline = topOutline(:,2) - min(topOutline(:,2));
o(2).v = [100 100]; o(2).xoutline = sideLeftOutline(:,1) - min(sideLeftOutline(:,1)); o(2).youtline = sideLeftOutline(:,2) - min(sideLeftOutline(:,2));
o(3).v = [100 100]; o(3).xoutline = sideRightOutline(:,1) - min(sideRightOutline(:,1)); o(3).youtline = sideRightOutline(:,2) - min(sideRightOutline(:,2));
nOutlines = length(o);

%% --- Generate Priority Optode Geometry ---
[detectors_left, sources_left, detectors_right, sources_right, ...
 det_centroids_left, src_centroids_left, det_centroids_right, src_centroids_right, ...
 all_sources_left, all_sources_right, all_detectors_left, all_detectors_right] = ...
 get_UHD_optode_pos_modified(sideLeftOutline, sideRightOutline);

conn_det_left = find_connections(det_centroids_left, sources_left, 3, 15);
conn_src_left = find_connections(src_centroids_left, sources_left, 3, 15);
all_connections_left_raw = [conn_det_left; conn_src_left];
conn_det_right = find_connections(det_centroids_right, sources_right, 3, 15);
conn_src_right = find_connections(src_centroids_right, sources_right, 3, 15);
all_connections_right_raw = [conn_det_right; conn_src_right];

exclusion_pts_left = [all_sources_left; all_detectors_left; det_centroids_left; src_centroids_left];
exclusion_pts_right = [all_sources_right; all_detectors_right; det_centroids_right; src_centroids_right];
offset_left = min(sideLeftOutline(:,1:2));
[v_conn_left, e_conn_left] = convert_plot_lines_to_graph(all_connections_left_raw);
if ~isempty(v_conn_left), v_conn_left = v_conn_left - offset_left; exclusion_pts_left = exclusion_pts_left - offset_left; all_connections_left_shifted = all_connections_left_raw; if ~isempty(all_connections_left_shifted), all_connections_left_shifted(~isnan(all_connections_left_shifted(:,1)),:) = all_connections_left_shifted(~isnan(all_connections_left_shifted(:,1)),:) - offset_left; end; end
offset_right = min(sideRightOutline(:,1:2));
[v_conn_right, e_conn_right] = convert_plot_lines_to_graph(all_connections_right_raw);
if ~isempty(v_conn_right), v_conn_right = v_conn_right - offset_right; exclusion_pts_right = exclusion_pts_right - offset_right; all_connections_right_shifted = all_connections_right_raw; if ~isempty(all_connections_right_shifted), all_connections_right_shifted(~isnan(all_connections_right_shifted(:,1)),:) = all_connections_right_shifted(~isnan(all_connections_right_shifted(:,1)),:) - offset_right; end; end

%% --- Set up Seams ---
iSeam = 1;
xx = o(iSeam).xoutline; yy = o(iSeam).youtline;
o(iSeam).seamLength = zeros(length(topOutline),1); o(iSeam).seamIs = zeros(length(topOutline),1); o(iSeam).seamIs(topIdx(1)) = 1;
for ii = topIdx(1)+1:topIdx(2), o(iSeam).seamLength(ii) = o(iSeam).seamLength(ii-1) + norm( [xx(ii)-xx(ii-1) yy(ii)-yy(ii-1)] ); o(iSeam).seamIs(ii) = 1; end
o(iSeam).seamIs(topIdx(3)) = 1;
for ii = topIdx(3)+1:topIdx(4), o(iSeam).seamLength(ii) = o(iSeam).seamLength(ii-1) + norm( [xx(ii)-xx(ii-1) yy(ii)-yy(ii-1)] ); o(iSeam).seamIs(ii) = 1; end
iSeam = 2;
xx = o(iSeam).xoutline; yy = o(iSeam).youtline;
o(iSeam).seamLength = zeros(length(sideLeftOutline),1); o(iSeam).seamIs = zeros(length(sideLeftOutline),1); o(iSeam).seamIs(sideLeftIdx(1)) = 1;
for ii = sideLeftIdx(1)+1 : sideLeftIdx(2), o(iSeam).seamLength(ii) = o(iSeam).seamLength(ii-1) + norm( [xx(ii)-xx(ii-1) yy(ii)-yy(ii-1)] ); o(iSeam).seamIs(ii) = 1; end
iSeam = 3;
xx = o(iSeam).xoutline; yy = o(iSeam).youtline;
o(iSeam).seamLength = zeros(length(sideRightOutline),1); o(iSeam).seamIs = zeros(length(sideRightOutline),1); o(iSeam).seamIs(sideRightIdx(1)) = 1;
for ii = sideRightIdx(1)+1 : sideRightIdx(2), o(iSeam).seamLength(ii) = o(iSeam).seamLength(ii-1) + norm( [xx(ii)-xx(ii-1) yy(ii)-yy(ii-1)] ); o(iSeam).seamIs(ii) = 1; end

%% --- Main Geometry Generation Loop ---
for iO = 1:nOutlines
    [xgrid,ygrid] = meshgrid(0:max(o(iO).xoutline)*1.05, 0:max(o(iO).youtline)*1.05);
    final_mask = inpolygon(xgrid,ygrid,o(iO).xoutline,o(iO).youtline);
    
    v_conn = []; e_conn = [];
    if iO == 2, final_mask = create_exclusion_mask(final_mask, xgrid, ygrid, exclusion_pts_left, all_connections_left_shifted, clearance); v_conn = v_conn_left; e_conn = e_conn_left;
    elseif iO == 3, final_mask = create_exclusion_mask(final_mask, xgrid, ygrid, exclusion_pts_right, all_connections_right_shifted, clearance); v_conn = v_conn_right; e_conn = e_conn_right;
    end
    
    [v_hex, e_hex, ~, ~] = fillCapWithHexagons_func(o(iO).v, hEdge, final_mask);
    
    num_hex_verts = size(v_hex, 1);
    temp_v = v_hex; temp_e = e_hex;
    if ~isempty(v_conn), e_conn_remapped = e_conn + num_hex_verts; temp_v = [v_hex; v_conn]; temp_e = [e_hex; e_conn_remapped]; end
    if iO > 1 && ~isempty(v_conn), stitching_edges = stitch_geometries(v_hex, e_hex, v_conn, num_hex_verts, connection_threshold); temp_e = [temp_e; stitching_edges]; end
    
    [v_clean, e_clean] = clean_mesh(temp_v, temp_e);
    o(iO).v = v_clean; o(iO).e = e_clean; o(iO).Imask = final_mask;

    outline_poly = [o(iO).xoutline, o(iO).youtline];
    [o(iO).vOut, o(iO).eOut] = recalculate_boundary_by_intersection(o(iO).v, o(iO).e, outline_poly);
end

%% DEBUG PLOT 1: Check initial geometry before final transformations
% Uncomment the block below to see the generated panels
%{
figure;
for iO = 1:nOutlines
    subplot(1, 3, iO);
    hold on;
    title(['Panel ' num2str(iO) ' (Before Final Stitching)']);
    
    % Plot the main mesh edges
    e_plot = o(iO).e;
    v_plot = o(iO).v;
    plot([v_plot(e_plot(:,1),1), v_plot(e_plot(:,2),1)]', [v_plot(e_plot(:,1),2), v_plot(e_plot(:,2),2)]', 'k-');
    
    % Plot the calculated outer boundary edges in red
    if isfield(o(iO), 'vOut') && ~isempty(o(iO).vOut)
        eOut_plot = o(iO).eOut;
        vOut_plot = o(iO).vOut;
        plot([vOut_plot(eOut_plot(:,1),1), vOut_plot(eOut_plot(:,2),1)]', [vOut_plot(eOut_plot(:,1),2), vOut_plot(eOut_plot(:,2),2)]', 'r-', 'LineWidth', 2);
    end
    
    % Plot the outline
    plot(o(iO).xoutline, o(iO).youtline, 'b--', 'LineWidth', 1.5);
    
    axis equal;
    hold off;
end
return; % Use return to stop execution here to check the plot
%}

%% --- All Downstream Processing (from your original script) ---
for iO = 2:3 % Reduced loop since iO=1 has no seam extensions
    if ~isfield(o(iO),'vOut') || isempty(o(iO).vOut), continue; end
    xx = o(iO).xoutline; yy = o(iO).youtline;
    xy1 = [xx yy [xx(2:end); xx(1)] [yy(2:end); yy(1)] ];
    
    eSeamX = zeros(size(o(iO).eOut,1),1); eSeamXidx = zeros(size(o(iO).eOut,1),1);
    eSeamXpts = zeros(size(o(iO).eOut,1),6); eSeamXlen = zeros(size(o(iO).eOut,1),1);
    eSeamXtheta = zeros(size(o(iO).eOut,1),1); eSeamXrot = zeros(size(o(iO).eOut,1),2,2);
    eSeamXptsRot = zeros(size(o(iO).eOut,1),6);
    
    [~, v_out_indices_in_v] = ismember(round(o(iO).vOut,4), round(o(iO).v,4), 'rows');

    for iE = 1:size(o(iO).eOut,1)
        v1_idx = v_out_indices_in_v(o(iO).eOut(iE,1));
        v2_idx = v_out_indices_in_v(o(iO).eOut(iE,2));
        if v1_idx == 0 || v2_idx == 0, continue; end
        xy2 = [o(iO).v(v1_idx,:) o(iO).v(v2_idx,:)];
        
        out = lineSegmentIntersect(xy1, xy2);
        idx = find(out.intAdjacencyMatrix==1);
        
        if ~isempty(idx)
            idx = idx(1);
            if o(iO).seamIs(idx)==1
                eSeamX(iE) = 1; eSeamXidx(iE) = idx;
                if inpolygon(o(iO).v(v1_idx,1), o(iO).v(v1_idx,2), xx, yy), inside_v_idx = v1_idx; else, inside_v_idx = v2_idx; end
                eSeamXpts(iE,1:4) = [o(iO).v(inside_v_idx,:) out.intMatrixX(idx) out.intMatrixY(idx)];
                
                safe_idx_end = min(idx+1, length(o(iO).seamLength));
                eSeamXlen(iE) = o(iO).seamLength(idx) + diff(o(iO).seamLength(idx:safe_idx_end)) * out.intNormalizedDistance1To2(idx);
                
                theta = atan2(xy1(idx,4)-xy1(idx,2), xy1(idx,3)-xy1(idx,1));
                eSeamXtheta(iE) = theta; R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
                if iO==2, pt = get_point_on_outline([(o(1).xoutline(1:topIdx(2))) (o(1).youtline(1:topIdx(2)))], eSeamXlen(iE));
                else, pt = get_point_on_outline([flipud(o(1).xoutline(topIdx(3):topIdx(4))) flipud(o(1).youtline(topIdx(3):topIdx(4)))], eSeamXlen(iE)); end
                eSeamXptsRot(iE,1:2) = (R * (eSeamXpts(iE,1:2)-eSeamXpts(iE,3:4))')' + pt;
                eSeamXptsRot(iE,3:4) = pt; eSeamXrot(iE,:,:) = R;
            end
        end
    end
    o(iO).eSeamX = eSeamX; o(iO).eSeamXidx = eSeamXidx; o(iO).eSeamXpts = eSeamXpts;
    o(iO).eSeamXlen = eSeamXlen; o(iO).eSeamXtheta = eSeamXtheta; o(iO).eSeamXrot = eSeamXrot; o(iO).eSeamXptsRot = eSeamXptsRot;
end

for iO = 2:3
    if isfield(o(iO), 'eSeamX') && any(o(iO).eSeamX)
        for iE = 1:length(o(iO).eSeamX)
            if o(iO).eSeamX(iE)==1
                p1 = o(iO).eSeamXptsRot(iE,3:4); p2lst = o(1).v;
                if isempty(p2lst), continue; end
                [~,idx] = min(sum(bsxfun(@minus, p2lst, p1).^2, 2));
                o(iO).eSeamXptsRot(iE,5:6) = p2lst(idx,:);
                p0 = p2lst(idx,:) - p1; theta = o(iO).eSeamXtheta(iE);
                R_inv = [cos(theta) sin(theta); -sin(theta) cos(theta)];
                o(iO).eSeamXpts(iE,5:6) = (R_inv * p0')' + o(iO).eSeamXpts(iE,3:4);
            end
        end
        lst = find(o(iO).eSeamX==1);
        o(iO).eSeamExt = o(iO).eSeamXpts(lst,[1 2 5 6]);
        o(iO).eSeamExtRot = o(iO).eSeamXptsRot(lst,[1 2 5 6]);
    end
end

for iO = 2:3
    if isfield(o(iO), 'eSeamExt') && ~isempty(o(iO).eSeamExt)
        o(iO).eExtensions = struct('pts',{});
        for pt = 1:size(o(iO).eSeamExt,1)
            [~, v_idx_in_v] = ismember(round(o(iO).eSeamExt(pt,1:2),4), round(o(iO).v,4), 'rows');
            if v_idx_in_v == 0, continue; end
            eidx = find(o(iO).e(:,1) == v_idx_in_v | o(iO).e(:,2) == v_idx_in_v);
            other_vidx = setdiff(o(iO).e(eidx,:), v_idx_in_v);
            pts = [];
            for u = 1:length(other_vidx), pts = [pts, o(iO).v(v_idx_in_v,:), o(iO).v(other_vidx(u),:)]; end
            o(iO).eExtensions(pt).pts = pts;
            Movingpts = [o(iO).eSeamExt(pt,1), o(iO).eSeamExt(pt,2); o(iO).eSeamExt(pt,3), o(iO).eSeamExt(pt,4)];
            FixedPts = [o(iO).eSeamExtRot(pt,1), o(iO).eSeamExtRot(pt,2); o(iO).eSeamExtRot(pt,3), o(iO).eSeamExtRot(pt,4)];
            if any(isnan(Movingpts(:))) || any(isnan(FixedPts(:))), continue; end
            tform = fitgeotrans(Movingpts, FixedPts, 'nonreflectivesimilarity');
            if isempty(pts), continue; end
            Xpts = pts(1:2:end); Ypts = pts(2:2:end);
            [Txpts,Typts] = transformPointsForward(tform, Xpts, Ypts);
            o(iO).eExtensionsRot(pt).pts = reshape([Txpts,Typts]', 1, []);
        end
    end
end

v_top = o(1).v; e_top = o(1).e;
for iO = 2:3
    if isfield(o(iO), 'eSeamExtRot') && ~isempty(o(iO).eSeamExtRot)
         o(iO).eExtensionsTop = struct('pts',{});
        for pt = 1:size(o(iO).eSeamExtRot,1)
            [~, v_idx_in_v] = ismember(round(o(iO).eSeamExtRot(pt,3:4),4), round(v_top,4), 'rows');
            if v_idx_in_v == 0, continue; end
            eidx = find(e_top(:,1) == v_idx_in_v | e_top(:,2) == v_idx_in_v);
            other_vidx = setdiff(e_top(eidx,:),v_idx_in_v);
            pts = [];
            for u = 1:length(other_vidx), pts = [pts, v_top(v_idx_in_v,:), v_top(other_vidx(u),:)]; end
            o(iO).eExtensionsTop(pt).pts = pts;
        end
    end
end

end

% =========================================================================
% LOCAL HELPER FUNCTIONS (Must be in the same .m file)
% =========================================================================

function [vOut, eOut] = recalculate_boundary_by_intersection(v_final, e_final, outline_poly)
% This robustly finds edges that physically intersect the outline.
    vOut = []; eOut = [];
    if isempty(v_final) || isempty(e_final), return; end
    xx = outline_poly(:,1); yy = outline_poly(:,2);
    outline_segments = [xx, yy, [xx(2:end); xx(1)], [yy(2:end); yy(1)]];
    crossing_edge_indices = [];
    for i = 1:size(e_final, 1)
        v1 = v_final(e_final(i, 1), :);
        v2 = v_final(e_final(i, 2), :);
        mesh_edge_segment = [v1, v2];
        out = lineSegmentIntersect(outline_segments, mesh_edge_segment);
        if any(out.intAdjacencyMatrix(:))
            crossing_edge_indices(end+1) = i;
        end
    end
    if isempty(crossing_edge_indices), return; end
    e_crossing = e_final(crossing_edge_indices, :);
    v_out_indices_in_v = unique(e_crossing(:));
    vOut = v_final(v_out_indices_in_v, :);
    [~, eOut] = ismember(e_crossing, v_out_indices_in_v);
end

function [v_new, e_new] = clean_mesh(v, e)
    if isempty(v) || isempty(e), v_new = v; e_new = e; return; end
    [~, I, J] = unique(round(v, 4), 'rows', 'first');
    map = zeros(size(v,1),1); map(I) = 1:length(I);
    e_new = map(J(e));
    v_new = v(I,:);
    e_new(any(e_new == 0, 2), :) = [];
    e_new = sort(e_new, 2);
    [~, Ia, ~] = unique(e_new, 'rows', 'first');
    e_new = e_new(Ia, :);
end

function stitching_edges = stitch_geometries(v_hex, e_hex, v_conn, num_hex_verts, threshold)
    stitching_edges = [];
    if isempty(v_hex) || isempty(v_conn), return; end
    edge_counts = histcounts(e_hex(:), 1:num_hex_verts+1);
    boundary_hex_indices = find(edge_counts > 0 & edge_counts < 6)'; % Accept any non-interior hex node
    if isempty(boundary_hex_indices), return; end
    v_hex_boundary = v_hex(boundary_hex_indices, :);
    [indices, dists] = knnsearch(v_hex_boundary, v_conn, 'K', 1);
    for i = 1:size(v_conn, 1)
        if dists(i) < threshold
            conn_v_idx = num_hex_verts + i;
            hex_v_idx = boundary_hex_indices(indices(i));
            stitching_edges = [stitching_edges; conn_v_idx, hex_v_idx];
        end
    end
end

function modified_mask = create_exclusion_mask(base_mask, xgrid, ygrid, ex_pts, ex_lines, clearance)
    if (isempty(ex_pts) || all(isnan(ex_pts(:)))) && (isempty(ex_lines) || all(isnan(ex_lines(:)))), modified_mask = base_mask; return; end
    grid_pts = [xgrid(base_mask), ygrid(base_mask)];
    if isempty(grid_pts), modified_mask = base_mask; return; end
    to_remove_mask = false(size(grid_pts, 1), 1);
    if ~isempty(ex_pts) && ~all(isnan(ex_pts(:)))
        min_dist_to_points = min(pdist2(grid_pts, ex_pts(~any(isnan(ex_pts),2),:)), [], 2);
        to_remove_mask = to_remove_mask | (min_dist_to_points < clearance);
    end
    if ~isempty(ex_lines) && ~all(isnan(ex_lines(:)))
        conn_segments = [];
        for i = 1:3:size(ex_lines, 1)
            segment = ex_lines(i:i+1, :);
            if ~any(isnan(segment), 'all'), conn_segments(end+1, :, :) = segment; end
        end
        if ~isempty(conn_segments)
            min_dist_to_lines = inf(size(grid_pts, 1), 1);
            for i = 1:size(conn_segments, 1)
                dists_to_segment = point_to_segment_distance_vectorized(grid_pts, squeeze(conn_segments(i, 1, :))', squeeze(conn_segments(i, 2, :))');
                min_dist_to_lines = min(min_dist_to_lines, dists_to_segment);
            end
            to_remove_mask = to_remove_mask | (min_dist_to_lines < clearance);
        end
    end
    modified_mask = base_mask;
    original_indices = find(base_mask);
    modified_mask(original_indices(to_remove_mask)) = 0;
end

function [v, e] = convert_plot_lines_to_graph(lines)
    if isempty(lines), v = []; e = []; return; end
    valid_lines = lines(~isnan(lines(:,1)), :);
    if isempty(valid_lines), v = []; e = []; return; end
    starts = valid_lines(1:2:end, :); ends = valid_lines(2:2:end, :);
    all_points = [starts; ends];
    [v, ~, ic] = unique(all_points, 'rows', 'stable');
    e = [ic(1:size(starts, 1)), ic(size(starts, 1)+1:end)];
end

function connections = find_connections(centers, positions, k, max_distance)
    connections = [];
    if isempty(centers) || isempty(positions), return; end
    for i = 1:size(centers, 1)
        distances = pdist2(centers(i, :), positions);
        [sorted_distances, sorted_indices] = sort(distances);
        closest_count = 0;
        for j = 1:length(sorted_distances)
            if sorted_distances(j) <= max_distance
                if closest_count < k, connections = [connections ; centers(i, :); positions(sorted_indices(j), :); [NaN NaN]]; closest_count = closest_count + 1; else, break; end
            else, break; end
        end
    end
end



function [outer_curve_x, outer_curve_y] = create_outer_curve(curve_x, curve_y, outer_curve_offset)
    dx = gradient(curve_x); dy = gradient(curve_y);
    normal_x = dy; normal_y = -dx;
    magnitude = sqrt(normal_x.^2 + normal_y.^2); magnitude(magnitude == 0) = 1;
    unit_normal_x = normal_x ./ magnitude; unit_normal_y = normal_y ./ magnitude;
    outer_curve_x = curve_x + unit_normal_x * outer_curve_offset;
    outer_curve_y = curve_y + unit_normal_y * outer_curve_offset;
end

function dists = point_to_segment_distance_vectorized(pts, v1, v2)
    v1v2 = v2 - v1; ptv1 = pts - v1; L2 = sum(v1v2.^2);
    if L2 == 0, dists = sqrt(sum(ptv1.^2, 2)); return; end
    t = (ptv1 * v1v2') / L2; t = max(0, min(1, t));
    projections = v1 + t .* v1v2;
    dists = sqrt(sum((pts - projections).^2, 2));
end

