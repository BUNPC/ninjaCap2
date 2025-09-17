function o = createHybridCapDesign(sideLeftOutline, sideLeftIdx, sideRightOutline, sideRightIdx, topOutline, topIdx)
% createHybridCapDesign - Fills cap outlines with a hybrid geometry.
%
% This function first places a specific optode and connection pattern on the
% side panels, then fills the remaining area with a hexagonal lattice.
% The top panel is filled entirely with a hexagonal lattice. Finally, all
% three panels are stitched together at the seams.

%% --- Initial Setup & Hexagon Parameters ---
o(1).v = [50 50]; % first vertex location
o(1).xoutline = topOutline(:,1) - min(topOutline(:,1));
o(1).youtline = topOutline(:,2) - min(topOutline(:,2));

o(2).v = [100 100]; % first vertex location
o(2).xoutline = sideLeftOutline(:,1) - min(sideLeftOutline(:,1));
o(2).youtline = sideLeftOutline(:,2) - min(sideLeftOutline(:,2));

o(3).v = [100 100]; % first vertex location
o(3).xoutline = sideRightOutline(:,1) - min(sideRightOutline(:,1));
o(3).youtline = sideRightOutline(:,2) - min(sideRightOutline(:,2));

hEdge = 10; % hexagon edge length
nOutlines = length(o);

%% --- Generate Priority Connections for Side Panels ---
% Get the specific geometry from the 'get_UHD_optode_pos' logic.
% Note: We only pass the seam portion of the outlines as specified.
[detectors_left, sources_left, detectors_right, sources_right, ...
 det_centroids_left, src_centroids_left, det_centroids_right, src_centroids_right, ~, ~] = ...
 get_UHD_optode_pos(sideLeftOutline(sideLeftIdx(1):sideLeftIdx(2),:), ...
                    sideRightOutline(sideRightIdx(1):sideRightIdx(2),:), ...
                    topOutline);

% Generate the connection lines for the LEFT side panel
conn_det_left = find_connections(det_centroids_left, sources_left, 3, 15);
conn_src_left = find_connections(src_centroids_left, sources_left, 3, 15);
all_connections_left_raw = [conn_det_left; conn_src_left];

% Generate the connection lines for the RIGHT side panel
conn_det_right = find_connections(det_centroids_right, sources_right, 3, 15);
conn_src_right = find_connections(src_centroids_right, sources_right, 3, 15);
all_connections_right_raw = [conn_det_right; conn_src_right];

% Convert connection lines to vertex/edge format and shift to local panel coordinates
offset_left = min(sideLeftOutline(:,1:2));
[v_conn_left, e_conn_left] = convert_plot_lines_to_graph(all_connections_left_raw);
v_conn_left = v_conn_left - offset_left;
all_connections_left_shifted = all_connections_left_raw;
all_connections_left_shifted(~isnan(all_connections_left_shifted)) = all_connections_left_shifted(~isnan(all_connections_left_shifted)) - offset_left(1);
all_connections_left_shifted(3:3:end, :) = NaN; % keep NaNs

offset_right = min(sideRightOutline(:,1:2));
[v_conn_right, e_conn_right] = convert_plot_lines_to_graph(all_connections_right_raw);
v_conn_right = v_conn_right - offset_right;
all_connections_right_shifted = all_connections_right_raw;
all_connections_right_shifted(~isnan(all_connections_right_shifted)) = all_connections_right_shifted(~isnan(all_connections_right_shifted)) - offset_right(1);
all_connections_right_shifted(3:3:end, :) = NaN; % keep NaNs

%% --- Set up Seams (Unchanged from original) ---
iSeam = 1;
xx = o(iSeam).xoutline;
yy = o(iSeam).youtline;
o(iSeam).seamLength = zeros(length(topOutline),1);
o(iSeam).seamIs = zeros(length(topOutline),1);
o(iSeam).seamIs(topIdx(1)) = 1;
for ii = topIdx(1)+1:topIdx(2)
    o(iSeam).seamLength(ii) = o(iSeam).seamLength(ii-1) + norm( [xx(ii)-xx(ii-1) yy(ii)-yy(ii-1)] );
    o(iSeam).seamIs(ii) = 1;
end
o(iSeam).seamIs(topIdx(3)) = 1;
for ii = topIdx(3)+1:topIdx(4)
    o(iSeam).seamLength(ii) = o(iSeam).seamLength(ii-1) + norm( [xx(ii)-xx(ii-1) yy(ii)-yy(ii-1)] );
    o(iSeam).seamIs(ii) = 1;
end

iSeam = 2;
o(iSeam).seamStartIdx = sideLeftIdx(1);
o(iSeam).seamEndIdx = sideLeftIdx(2);
o(iSeam).seamLength = zeros(length(sideLeftOutline),1);
o(iSeam).seamIs = zeros(length(sideLeftOutline),1);
xx = o(iSeam).xoutline;
yy = o(iSeam).youtline;
o(iSeam).seamIs(o(iSeam).seamStartIdx) = 1;
for ii = o(iSeam).seamStartIdx+1 : o(iSeam).seamEndIdx
    o(iSeam).seamLength(ii) = o(iSeam).seamLength(ii-1) + norm( [xx(ii)-xx(ii-1) yy(ii)-yy(ii-1)] );
    o(iSeam).seamIs(ii) = 1;
end

iSeam = 3;
o(iSeam).seamStartIdx = sideRightIdx(1);
o(iSeam).seamEndIdx = sideRightIdx(2);
o(iSeam).seamLength = zeros(length(sideRightOutline),1);
o(iSeam).seamIs = zeros(length(sideRightOutline),1);
xx = o(iSeam).xoutline;
yy = o(iSeam).youtline;
o(iSeam).seamIs(o(iSeam).seamStartIdx) = 1;
for ii = o(iSeam).seamStartIdx+1 : o(iSeam).seamEndIdx
    o(iSeam).seamLength(ii) = o(iSeam).seamLength(ii-1) + norm( [xx(ii)-xx(ii-1) yy(ii)-yy(ii-1)] );
    o(iSeam).seamIs(ii) = 1;
end

%% --- Create Mask for Each Panel ---
for iO = 1:nOutlines
    [xgrid,ygrid] = meshgrid(1:max(o(iO).xoutline)*1.05,1:max(o(iO).youtline)*1.05);
    o(iO).Imask = inpolygon(xgrid,ygrid,o(iO).xoutline,o(iO).youtline);
end

%% --- MODIFIED SECTION: Fill Mask with Hexagons and Merge ---
for iO = 1:nOutlines
    % First, generate a full hexagonal grid for the panel's shape
    % We assume a function 'fillCapWithHexagons_func' exists as in the original code.
    % If you don't have it, this would be a function that generates a hex grid inside a mask.
    [v_hex, e_hex, vOut_hex, eOut_hex] = fillCapWithHexagons_func( o(iO).v, hEdge, o(iO).Imask );
    
    if iO == 1 % Top Panel: Use only the hexagonal fill
        o(iO).v = v_hex;
        o(iO).e = e_hex;
        o(iO).vOut = vOut_hex;
        o(iO).eOut = eOut_hex;
        
    else % Side Panels: Prune hex grid and merge with connections
        if iO == 2 % Left Panel
            v_conn = v_conn_left;
            e_conn = e_conn_left;
            conn_lines_shifted = all_connections_left_shifted;
        else % iO == 3, Right Panel
            v_conn = v_conn_right;
            e_conn = e_conn_right;
            conn_lines_shifted = all_connections_right_shifted;
        end
        
        % Prune hex edges that are too close to the priority connection struts
        e_hex_pruned = prune_hex_edges(v_hex, e_hex, conn_lines_shifted, hEdge / 2);
        
        % Combine the pruned hex grid with the connection grid
        num_hex_verts = size(v_hex, 1);
        e_conn_remapped = e_conn + num_hex_verts;
        
        o(iO).v = [v_hex; v_conn];
        o(iO).e = [e_hex_pruned; e_conn_remapped];
        
        % The outer boundary for stitching is still determined by the original hex fill
        o(iO).vOut = vOut_hex;
        o(iO).eOut = eOut_hex;
    end
end

%% --- Stitching and Cleanup (Largely unchanged from original) ---

% For outline edge (not outline seam), pull outside vertices to the outline
for iO = 1:length(o)
    xx = o(iO).xoutline;
    yy = o(iO).youtline;
    xy1 = [xx yy [xx(2:end); xx(1)] [yy(2:end); yy(1)] ];
    
    for iE = 1:size(o(iO).eOut,1)
        xy2 = [o(iO).vOut(o(iO).eOut(iE,1),:) o(iO).vOut(o(iO).eOut(iE,2),:)];
        out = lineSegmentIntersect(xy1, xy2); % Assumes lineSegmentIntersect.m is available
        idx = find(out.intAdjacencyMatrix==1);
        if ~isempty(idx) && o(iO).seamIs(idx(1))==0
            o(iO).v(end+1,:) = o(iO).vOut(o(iO).eOut(iE,1),:);
            o(iO).v(end+1,:) = [out.intMatrixX(idx(1)) out.intMatrixY(idx(1))];
            o(iO).e(end+1,:) = [size(o(iO).v,1)-1 size(o(iO).v,1)];
        end
    end
end

% ... (The rest of your original code for stitching and cleanup follows) ...
% The sections "Work on struts crossing the seam", "move side panel
% extended edges", "Remove duplicates nodes", and the final loop for
% "eExtensions" are now applied to the hybrid geometry stored in 'o'.
% The "Remove duplicates" step is especially important for merging vertices.

%% Work on struts crossing the seam
for iO = 2:3
    v = o(iO).v; vOut = o(iO).vOut; e = o(iO).e; eOut = o(iO).eOut;
    xx = o(iO).xoutline; yy = o(iO).youtline;
    xy1 = [xx yy [xx(2:end); xx(1)] [yy(2:end); yy(1)] ];
    
    eSeamX = zeros(size(eOut,1),1); eSeamXidx = zeros(size(eOut,1),1);
    eSeamXpts = zeros(size(eOut,1),6); eSeamXlen = zeros(size(eOut,1),1);
    eSeamXtheta = zeros(size(eOut,1),1); eSeamXrot = zeros(size(eOut,1),2,2);
    eSeamXptsRot = zeros(size(eOut,1),6);
    
    for iE = 1:size(eOut,1)
        idx = []; foo=0;
        while isempty(idx)
            xy2 = [vOut(eOut(iE,1),:)+foo*(vOut(eOut(iE,1),:)-vOut(eOut(iE,2),:)) vOut(eOut(iE,2),:)+foo*(vOut(eOut(iE,2),:)-vOut(eOut(iE,1),:))];
            foo = foo + 0.1;
            out = lineSegmentIntersect(xy1, xy2);
            idx = find(out.intAdjacencyMatrix==1);
        end
        idx = idx(1);
        eSeamXidx(iE) = idx;
        if o(iO).seamIs(idx)==1
            eSeamX(iE) = 1;
            eSeamXpts(iE,1:4) = [vOut(eOut(iE,1),:) out.intMatrixX(idx) out.intMatrixY(idx)];
            eSeamXlen(iE) = o(iO).seamLength(idx)*(1-out.intNormalizedDistance1To2(idx)) + o(iO).seamLength(idx+1)*out.intNormalizedDistance1To2(idx);
            
            theta = asin( (xy1(idx,3)-xy1(idx,1)) / norm(xy1(idx,3:4) - xy1(idx,1:2)) );
            eSeamXtheta(iE) = theta;
            R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
            if iO==2 % left side
                pt = get_point_on_outline([(o(1).xoutline(1:topIdx(2))) (o(1).youtline(1:topIdx(2)))], eSeamXlen(iE));
            else % right side
                pt = get_point_on_outline([flipud(o(1).xoutline(topIdx(3):topIdx(4))) flipud(o(1).youtline(topIdx(3):topIdx(4)))], eSeamXlen(iE));
            end
            eSeamXptsRot(iE,1:2) = R * (eSeamXpts(iE,1:2)-eSeamXpts(iE,3:4))' + pt';
            eSeamXptsRot(iE,3:4) = pt';
            eSeamXrot(iE,:,:) = R;
        end
    end
    o(iO).eSeamX = eSeamX; o(iO).eSeamXidx = eSeamXidx; o(iO).eSeamXpts = eSeamXpts;
    o(iO).eSeamXlen = eSeamXlen; o(iO).eSeamXtheta = eSeamXtheta; o(iO).eSeamXrot = eSeamXrot;
    o(iO).eSeamXptsRot = eSeamXptsRot;
end

%% move side panel extended edges to nearest top panel vertex
for iO = 2:3
    for iE = 1:length(o(iO).eSeamX)
        if o(iO).eSeamX(iE)==1
            p1 = o(iO).eSeamXptsRot(iE,[3:4]);
            p2lst = o(1).vOut(o(1).eOut(:,1),:);
            rho = sum((ones(size(p2lst,1),1)*p1 - p2lst).^2,2).^0.5;
            [~,idx] = min(rho);
            
            o(iO).eSeamXptsRot(iE,[5:6]) = p2lst(idx,:);
            
            p0 = p2lst(idx,:) - p1;
            theta = o(iO).eSeamXtheta(iE);
            R = [cos(theta) sin(theta); -sin(theta) cos(theta)];
            o(iO).eSeamXpts(iE,[5:6]) = R * p0' + o(iO).eSeamXpts(iE,3:4)';
        end
    end
    lst = find(o(iO).eSeamX==1);
    o(iO).eSeamExt = o(iO).eSeamXpts(lst,[1 2 5 6]);
    o(iO).eSeamExtRot = o(iO).eSeamXptsRot(lst,[1 2 5 6]);
end

%% Remove duplicates nodes in v, e, vout and eout
for iO = 1:3
    v = o(iO).v;
    e = o(iO).e;
    [~,I,J] = unique(round(v, 4), 'rows', 'first'); % Round to handle floating point inaccuracies
    
    map = zeros(size(v,1),1);
    map(I) = 1:length(I);
    e_new = map(J(e));
    v_new = v(I,:);
    
    e_new(any(e_new == 0, 2), :) = []; % remove edges with zero index
    
    e_new = sort(e_new, 2);
    [~, Ia, ~] = unique(e_new, 'rows', 'first');
    e_new = e_new(Ia, :);
    
    o(iO).v = v_new;
    o(iO).e = e_new;
end


%% Final eExtensions sections (Unchanged)
for iO = 2:3
    v = o(iO).v; e = o(iO).e;
    for pt = 1:size(o(iO).eSeamExt,1)
        vidx = find(v(:,1) == o(iO).eSeamExt(pt,1) & v(:,2) == o(iO).eSeamExt(pt,2));
        if isempty(vidx), continue, end
        eidx = find(e(:,1) == vidx | e(:,2) == vidx);
        
        pts = [];
        vidxs = e(eidx,:);
        other_vidx = setdiff(vidxs(:),vidx);
        for u =1:length(other_vidx)
            pts = [pts v(vidx,:) v(other_vidx(u),:)];
        end
        o(iO).eExtensions(pt).pts = pts;

        Movingpts = [o(iO).eSeamExt(pt,1) o(iO).eSeamExt(pt,2); o(iO).eSeamExt(pt,3) o(iO).eSeamExt(pt,4)];
        FixedPts = [o(iO).eSeamExtRot(pt,1) o(iO).eSeamExtRot(pt,2); o(iO).eSeamExtRot(pt,3) o(iO).eSeamExtRot(pt,4)];
        tform = fitgeotrans(Movingpts,FixedPts,'nonreflectivesimilarity');
        Xpts = o(iO).eExtensions(pt).pts(1:2:end);
        Ypts = o(iO).eExtensions(pt).pts(2:2:end);
        [Txpts,Typts] = transformPointsForward(tform, Xpts, Ypts);
        temppts = [Txpts;Typts];
        o(iO).eExtensionsRot(pt).pts = temppts(:)';
    end
end

v = o(1).v; e = o(1).e;
for iO = 2:3
    for pt = 1:size(o(iO).eSeamExtRot,1)
        vidx = find(v(:,1) == o(iO).eSeamExtRot(pt,3) & v(:,2) == o(iO).eSeamExtRot(pt,4));
        if isempty(vidx), continue, end
        eidx = find(e(:,1) == vidx | e(:,2) == vidx);
        
        pts = [];
        vidxs = e(eidx,:);
        other_vidx = setdiff(vidxs(:),vidx);
        for u =1:length(other_vidx)
            pts = [pts v(vidx,:) v(other_vidx(u),:)];
        end
        o(iO).eExtensionsTop(pt).pts  = pts;
        
        FixedPts = [o(iO).eSeamExt(pt,1) o(iO).eSeamExt(pt,2); o(iO).eSeamExt(pt,3) o(iO).eSeamExt(pt,4)];
        Movingpts = [o(iO).eSeamExtRot(pt,1) o(iO).eSeamExtRot(pt,2); o(iO).eSeamExtRot(pt,3) o(iO).eSeamExtRot(pt,4)];
        tform = fitgeotrans(Movingpts,FixedPts,'nonreflectivesimilarity');
        Xpts = o(iO).eExtensionsTop(pt).pts(1:2:end);
        Ypts = o(iO).eExtensionsTop(pt).pts(2:2:end);
        [Txpts,Typts] = transformPointsForward(tform, Xpts, Ypts);
        temppts = [Txpts;Typts];
        o(iO).eExtensionsTopRot(pt).pts = temppts(:)';
    end
end

end % End of main function

% -------------------------------------------------------------------------
% LOCAL HELPER FUNCTIONS
% (Place these in the same .m file, below the main function)
% -------------------------------------------------------------------------

function e_pruned = prune_hex_edges(v_hex, e_hex, conn_lines, threshold)
% Removes edges from e_hex if their midpoint is too close to any conn_line
    
    % Extract individual line segments from the connection plot data
    conn_segments = [];
    for i = 1:3:size(conn_lines, 1)
        segment = conn_lines(i:i+1, :);
        if ~any(isnan(segment(:)))
            conn_segments(end+1, :, :) = segment;
        end
    end

    if isempty(conn_segments)
        e_pruned = e_hex;
        return;
    end
    
    to_keep = true(size(e_hex, 1), 1);
    for i = 1:size(e_hex, 1)
        p1 = v_hex(e_hex(i, 1), :);
        p2 = v_hex(e_hex(i, 2), :);
        midpoint = (p1 + p2) / 2;
        
        min_dist_to_any_segment = inf;
        for j = 1:size(conn_segments, 1)
            v1 = squeeze(conn_segments(j, 1, :))';
            v2 = squeeze(conn_segments(j, 2, :))';
            dist = point_to_segment_distance(midpoint, v1, v2);
            if dist < min_dist_to_any_segment
                min_dist_to_any_segment = dist;
            end
        end
        
        if min_dist_to_any_segment < threshold
            to_keep(i) = false;
        end
    end
    
    e_pruned = e_hex(to_keep, :);
end

function [v, e] = convert_plot_lines_to_graph(lines)
% Converts a list of line segments for plotting into a vertex and edge list
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

% --- Include the helper functions you provided, like get_UHD_optode_pos,
% --- find_connections, point_to_segment_distance, etc., here as local functions.

function [detectors_left_side, sources_left_side, detectors_right_side, ...
    sources_right_side, det_centroids_left, src_centroids_left, ...
    det_centroids_right, src_centroids_right, outer_sources_left_side, outer_sources_right_side] = ...
    get_UHD_optode_pos(sideLeftOutline, sideRightOutline, topOutline)
    % This is the get_UHD_optode_pos function you provided.
    % Pasted here for completeness.
    %% ---Fill optodes on side Left ---
    vertical_spacing = 12;
    outer_curve_factor = 1;
    outer_curve_offset = outer_curve_factor*vertical_spacing;
    for u = 1:2
        if  u == 1
        % Horizontal spacing between columns
            horizontal_spacing = vertical_spacing * sind(60);
            curve_y = sideLeftOutline(:,2); 
            curve_x = sideLeftOutline(:,1);
            [outer_curve_x, outer_curve_y] = create_outer_curve(curve_x, curve_y, outer_curve_offset);
        else
            horizontal_spacing = -vertical_spacing * sind(60);
            curve_y = sideRightOutline(:,2); 
            curve_x = sideRightOutline(:,1);
            [outer_curve_x, outer_curve_y] = create_outer_curve(curve_x, curve_y, -outer_curve_offset);
        end

        detectors = []; sources = [];
        detectors_all = []; sources_all = [];

        polygon_x = [curve_x(1); curve_x; curve_x(end)];
        polygon_y = [curve_y(1); curve_y; curve_y(end)];
        outer_polygon_x = [outer_curve_x(1); outer_curve_x; outer_curve_x(end)];
        outer_polygon_y = [outer_curve_y(1); outer_curve_y; outer_curve_y(end)];

        x_pos = curve_x(end); 
        y_start_col1 = []; 
        for y_pos = max(curve_y)+outer_curve_factor*vertical_spacing:-vertical_spacing:min(curve_y)-outer_curve_factor*vertical_spacing
            detectors = [detectors; x_pos, y_pos];
            if isempty(y_start_col1)
               y_start_col1 = y_pos;
            end
            detectors_all = [detectors_all; x_pos, y_pos];
        end
        detectors = detectors(outer_curve_factor+1:end-outer_curve_factor,:);

        column_is_detector_only = false; 
        x_pos = x_pos + horizontal_spacing;
        if u == 1
            condition = @(x_pos) x_pos < max(outer_curve_x);
        else
            condition = @(x_pos) x_pos > min(outer_curve_x);
        end
        det_or_src = 0;
        while condition(x_pos)
            if column_is_detector_only
                for y_pos = max(curve_y)+outer_curve_factor*vertical_spacing:-vertical_spacing:min(curve_y)-outer_curve_factor*vertical_spacing
                    if inpolygon(x_pos, y_pos, polygon_x, polygon_y)
                        detectors = [detectors; x_pos, y_pos];
                    end
                    if inpolygon(x_pos, y_pos, outer_polygon_x, outer_polygon_y) || inpolygon(x_pos, y_pos, polygon_x, polygon_y)
                       detectors_all = [detectors_all; x_pos, y_pos];
                    end
                end
            else
                y_start_col2 = y_start_col1 +(outer_curve_factor-0.5)*vertical_spacing;
                for i = 0:25
                    if mod(outer_curve_factor+det_or_src,2) == 1
                        detector_y = y_start_col2 - i * (2 * vertical_spacing);
                        source_y = detector_y - vertical_spacing;
                    else
                        source_y = y_start_col2 - i * (2 * vertical_spacing);
                        detector_y = source_y - vertical_spacing;
                    end
                    if inpolygon(x_pos, detector_y, polygon_x, polygon_y)
                        detectors = [detectors; x_pos, detector_y]; 
                    end
                    if inpolygon(x_pos, source_y, polygon_x, polygon_y)
                        sources = [sources; x_pos, source_y]; 
                    end
                    if inpolygon(x_pos, detector_y, outer_polygon_x, outer_polygon_y) || inpolygon(x_pos, detector_y, polygon_x, polygon_y)
                        detectors_all = [detectors_all; x_pos, detector_y];
                    end
                    if inpolygon(x_pos, source_y, outer_polygon_x, outer_polygon_y) || inpolygon(x_pos, source_y, polygon_x, polygon_y)
                        sources_all = [sources_all; x_pos, source_y];
                    end
                end
                det_or_src = det_or_src+1;
            end
            x_pos = x_pos + horizontal_spacing;
            column_is_detector_only = ~column_is_detector_only;
        end
        if  u == 1
            detectors_left_side = detectors;
            sources_left_side = sources;
            outer_detectors_left_side = detectors_all;
            outer_sources_left_side = sources_all;
        else
            detectors_right_side = detectors;
            sources_right_side = sources;
            outer_detectors_right_side = detectors_all;
            outer_sources_right_side = sources_all;
        end
    end

    det_centroids_left = get_det_centroids(outer_detectors_left_side,12, 15);
    det_centroids_right = get_det_centroids(outer_detectors_right_side,12, 15);
    src_centroids_left = det_centroids_left+[12*sind(60) 0];
    src_centroids_right = det_centroids_right-[12*sind(60) 0];
end

function connections = find_connections(centers, positions, k, max_distance)
% This is the find_connections function you provided.
    connections = []; 
    for i = 1:size(centers, 1)
        center = centers(i, :);
        distances = pdist2(center, positions);
        [sorted_distances, sorted_indices] = sort(distances);
        closest_count = 0;
        for j = 1:length(sorted_distances)
            if sorted_distances(j) <= max_distance
                if closest_count < k
                    position_index = sorted_indices(j);
                    connections = [connections ; center; positions(position_index, :); [NaN NaN]];
                    closest_count = closest_count + 1;
                else
                    break;
                end
            else
                break;
            end
        end
    end
end

function dist = point_to_segment_distance(pt, v1, v2)
% This is the point_to_segment_distance function you provided.
    v1v2 = v2 - v1;
    ptv1 = pt - v1;
    L2 = sum(v1v2.^2);
    if L2 == 0
        dist = norm(ptv1);
        return;
    end
    t = dot(ptv1, v1v2) / L2;
    t = max(0, min(1, t));
    projection = v1 + t * v1v2;
    dist = norm(pt - projection);
end

function [outer_curve_x, outer_curve_y] = create_outer_curve(curve_x, curve_y, outer_curve_offset)
% This is the create_outer_curve function you provided.
    dx = gradient(curve_x);
    dy = gradient(curve_y);
    normal_x = dy;
    normal_y = -dx;
    magnitude = sqrt(normal_x.^2 + normal_y.^2);
    magnitude(magnitude == 0) = 1; 
    unit_normal_x = normal_x ./ magnitude;
    unit_normal_y = normal_y ./ magnitude;
    outer_curve_x = curve_x + unit_normal_x * outer_curve_offset;
    outer_curve_y = curve_y + unit_normal_y * outer_curve_offset;
end

% function det_centroids = get_det_centroids(detectors, spacing, tolerance)
% % Dummy function for get_det_centroids since it wasn't provided.
% % This one finds equilateral triangles of detectors pointing up or down.
%     det_centroids = [];
%     tol_sq = tolerance^2;
%     
%     % Use pdist2 for efficient distance calculations
%     dists = squareform(pdist(detectors));
%     
%     for i = 1:size(detectors,1)
%         % Find neighbors at the correct spacing
%         neighbors_idx = find(abs(dists(i,:) - spacing) < spacing*0.1);
%         
%         if length(neighbors_idx) >= 2
%             % Check pairs of neighbors
%             combos = nchoosek(neighbors_idx, 2);
%             for j = 1:size(combos,1)
%                 p2_idx = combos(j,1);
%                 p3_idx = combos(j,2);
%                 
%                 % Check if these three points form an equilateral triangle
%                 dist_23 = norm(detectors(p2_idx,:) - detectors(p3_idx,:));
%                 
%                 if abs(dist_23 - spacing) < spacing*0.1
%                     centroid = mean(detectors([i, p2_idx, p3_idx],:));
%                     det_centroids = [det_centroids; centroid];
%                 end
%             end
%         end
%     end
%     % Remove duplicate centroids
%     if ~isempty(det_centroids)
%         det_centroids = unique(round(det_centroids, 2), 'rows');
%     end
% end
% NOTE: You will also need to have your 'fillCapWithHexagons_func.m',
% 'lineSegmentIntersect.m', and 'get_point_on_outline.m' functions available
% on your MATLAB path for this to run. I have created a placeholder for
% 'get_det_centroids' as it was not included in your prompt.