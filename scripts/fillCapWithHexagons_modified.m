function o = fillCapWithHexagons_modified(sideLeftOutline, sideLeftIdx, sideRightOutline, sideRightIdx, topOutline, topIdx)
% This is the main function.

% --- Initial Setup ---
o(1).v = [50 50]; 
o(1).xoutline = topOutline(:,1) - min(topOutline(:,1));
o(1).youtline = topOutline(:,2) - min(topOutline(:,2));

o(2).v = [100 100];
o(2).xoutline = sideLeftOutline(:,1) - min(sideLeftOutline(:,1));
o(2).youtline = sideLeftOutline(:,2) - min(sideLeftOutline(:,2));

o(3).v = [100 100];
o(3).xoutline = sideRightOutline(:,1) - min(sideRightOutline(:,1));
o(3).youtline = sideRightOutline(:,2) - min(sideRightOutline(:,2));

hEdge = 10; 
nOutlines = length(o);


%% --- Generate Pre-defined Connections for Side Panels ---
[~, sources_left_side, ~, sources_right_side, det_centroids_left, src_centroids_left, det_centroids_right, src_centroids_right, ~, ~] = ...
    get_UHD_optode_pos(sideLeftOutline(sideLeftIdx(1):sideLeftIdx(2),:), ...
                       sideRightOutline(sideRightIdx(1):sideRightIdx(2),:), ...
                       topOutline);

% Generate connection lists
connections_left_det = find_connections(det_centroids_left, sources_left_side, 3, 15);
connections_left_src = find_connections(src_centroids_left, sources_left_side, 3, 15);
all_connections_left = [connections_left_det; connections_left_src];

connections_right_det = find_connections(det_centroids_right, sources_right_side, 3, 15);
connections_right_src = find_connections(src_centroids_right, sources_right_side, 3, 15);
all_connections_right = [connections_right_det; connections_right_src];

% Normalize connection coordinates
min_left = min(sideLeftOutline(:,:), [], 1);
all_connections_left(~any(isnan(all_connections_left),2),:) = all_connections_left(~any(isnan(all_connections_left),2),:) - min_left;

min_right = min(sideRightOutline(:,:), [], 1);
all_connections_right(~any(isnan(all_connections_right),2),:) = all_connections_right(~any(isnan(all_connections_right),2),:) - min_right;

% Convert to V/E format
[v_conn_left, e_conn_left] = connectionsToVE(all_connections_left);
[v_conn_right, e_conn_right] = connectionsToVE(all_connections_right);


%% --- Set up Seams ---
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


%% --- Fill Each Panel ---
for iO = 1:nOutlines
    
    % 1. Create a binary mask for the ENTIRE panel
    max_x = ceil(max(o(iO).xoutline) * 1.05);
    max_y = ceil(max(o(iO).youtline) * 1.05);
    if max_x < 1, max_x = 1; end
    if max_y < 1, max_y = 1; end
    [xgrid, ygrid] = meshgrid(1:max_x, 1:max_y);
    panel_mask = inpolygon(xgrid, ygrid, o(iO).xoutline, o(iO).youtline);
    
    % 2. Fill the ENTIRE panel mask with a hexagonal grid
    [v_hex, e_hex, vOut_hex, eOut_hex] = fillCapWithHexagons_func(o(iO).v, hEdge, panel_mask);
    
    % 3. Get the connection struts for this panel (if any)
    v_conn = []; 
    e_conn = [];
    connections_list = [];
    if iO == 2 && ~isempty(v_conn_left)
        v_conn = v_conn_left;
        e_conn = e_conn_left;
        connections_list = all_connections_left;
    elseif iO == 3 && ~isempty(v_conn_right)
        v_conn = v_conn_right;
        e_conn = e_conn_right;
        connections_list = all_connections_right;
    end
    
    % 4. Prune the hexagon grid to make space for the connections
    if ~isempty(connections_list)
        % Calculate the distance from each hex vertex to the nearest connection strut
        dist_to_connections = inf(size(v_hex, 1), 1);
        for k = 1:3:size(connections_list, 1)
            p1 = connections_list(k, :);
            p2 = connections_list(k+1, :);
            if any(isnan(p1)) || any(isnan(p2)), continue; end
            
            % Vectorized distance calculation for all hex vertices to this one segment
            v1v2 = p2 - p1;
            ptv1 = v_hex - p1;
            L2 = sum(v1v2.^2);
            if L2 == 0
                dists_to_segment = sqrt(sum(ptv1.^2, 2));
            else
                t = dot(ptv1, repmat(v1v2, size(ptv1,1), 1), 2) / L2;
                t = max(0, min(1, t));
                projection = p1 + t .* v1v2;
                dists_to_segment = sqrt(sum((v_hex - projection).^2, 2));
            end
            dist_to_connections = min(dist_to_connections, dists_to_segment);
        end
        
        % Identify hex vertices to KEEP (those far enough away from connections)
        keep_threshold = hEdge * 0.5; % This is the keep-out radius around connections
        hex_verts_to_keep_idx = find(dist_to_connections > keep_threshold);
        
        % Create a map from old hex vertex indices to new (pruned) indices
        hex_map = zeros(size(v_hex, 1), 1);
        hex_map(hex_verts_to_keep_idx) = 1:length(hex_verts_to_keep_idx);
        
        % Filter the hex vertices and edges
        v_hex_pruned = v_hex(hex_verts_to_keep_idx, :);
        e_hex_mapped = hex_map(e_hex);
        e_hex_pruned = e_hex_mapped(all(e_hex_mapped > 0, 2), :);
        
        % Overwrite the original hex grid with the pruned version
        v_hex = v_hex_pruned;
        e_hex = e_hex_pruned;
    end
    
    % 5. Combine the connection struts and the pruned hexagon grid
    num_conn_verts = size(v_conn, 1);
    v_combined = [v_conn; v_hex];
    e_combined = [e_conn; e_hex + num_conn_verts];
    
    if ~isempty(v_conn) && ~isempty(v_hex)
        [indices, dists] = dsearchn(v_hex, v_conn);
        connection_edge_threshold = hEdge * 1.7;
        
        new_edges = [];
        for i = 1:length(indices)
            if dists(i) < connection_edge_threshold
                new_edges = [new_edges; i, num_conn_verts + indices(i)];
            end
        end
        e_combined = [e_combined; new_edges];
    end
    
    % 6. Assign the final structure to the output object
    o(iO).v = v_combined;
    o(iO).e = e_combined;
    o(iO).vOut = vOut_hex;
    o(iO).eOut = eOut_hex;
end


%% --- Post-Processing ---
% Pull outside vertices to the outline
for iO = 1:length(o)
    if isempty(o(iO).vOut) || isempty(o(iO).eOut), continue; end
    xx = o(iO).xoutline;
    yy = o(iO).youtline;
    xy1 = [xx yy [xx(2:end); xx(1)] [yy(2:end); yy(1)] ];
    for iE = 1:size(o(iO).eOut,1)
        xy2 = [o(iO).vOut(o(iO).eOut(iE,1),:) o(iO).vOut(o(iO).eOut(iE,2),:)];
        out = lineSegmentIntersect(xy1, xy2);
        idx = find(out.intAdjacencyMatrix==1);
        if ~isempty(idx) && o(iO).seamIs(idx(1))==0
            o(iO).v(end+1,:) = o(iO).vOut(o(iO).eOut(iE,1),:);
            o(iO).v(end+1,:) = [out.intMatrixX(idx(1)) out.intMatrixY(idx(1))];
            o(iO).e(end+1,:) = [size(o(iO).v,1)-1 size(o(iO).v,1)];
        end
    end
end

%% --- Process Seam-Crossing Struts (Combined Loop) ---
for iO = 2:3
    if isempty(o(iO).vOut) || isempty(o(iO).eOut), continue; end
    
    % Part 1: Find intersections with seam
    vOut = o(iO).vOut;
    eOut = o(iO).eOut;
    xx = o(iO).xoutline;
    yy = o(iO).youtline;
    xy1 = [xx yy [xx(2:end); xx(1)] [yy(2:end); yy(1)] ];
    
    eSeamX = zeros(size(eOut,1),1);
    eSeamXidx = zeros(size(eOut,1),1);
    eSeamXpts = zeros(size(eOut,1),6);
    eSeamXlen = zeros(size(eOut,1),1);
    eSeamXtheta = zeros(size(eOut,1),1);
    eSeamXrot = zeros(size(eOut,1),2,2);
    eSeamXptsRot = zeros(size(eOut,1),6);
    
    for iE = 1:size(eOut,1)
        idx = []; foo=0;
        while isempty(idx) && foo < 10 
            xy2 = [vOut(eOut(iE,1),:)+foo*(vOut(eOut(iE,1),:)-vOut(eOut(iE,2),:)) vOut(eOut(iE,2),:)+foo*(vOut(eOut(iE,2),:)-vOut(eOut(iE,1),:))];
            foo = foo + 0.1;
            out = lineSegmentIntersect(xy1, xy2);
            idx = find(out.intAdjacencyMatrix==1);
        end
        if isempty(idx), continue; end 
        idx = idx(1);
        eSeamXidx(iE) = idx;
        if o(iO).seamIs(idx)==1
            eSeamX(iE) = 1;
            eSeamXpts(iE,1:4) = [vOut(eOut(iE,1),:) out.intMatrixX(idx) out.intMatrixY(idx)];
            eSeamXlen(iE) = o(iO).seamLength(idx)*(1-out.intNormalizedDistance1To2(idx)) + o(iO).seamLength(idx+1)*out.intNormalizedDistance1To2(idx);
            theta = asin( (xy1(idx,3)-xy1(idx,1)) / (norm(xy1(idx,3:4) - xy1(idx,1:2)) + eps) );
            eSeamXtheta(iE) = theta;
            if iO==2 % left side
                R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
                pt = get_point_on_outline([(o(1).xoutline(1:topIdx(2))) (o(1).youtline(1:topIdx(2)))], eSeamXlen(iE));
                eSeamXptsRot(iE,1:2) = R * (eSeamXpts(iE,1:2)-eSeamXpts(iE,3:4))' + pt';
                eSeamXptsRot(iE,3:4) = pt';
            else % right side
                R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
                pt = get_point_on_outline([flipud(o(1).xoutline(topIdx(3):topIdx(4))) flipud(o(1).youtline(topIdx(3):topIdx(4)))], eSeamXlen(iE));
                eSeamXptsRot(iE,1:2) = R * (eSeamXpts(iE,1:2)-eSeamXpts(iE,3:4))' + pt';
                eSeamXptsRot(iE,3:4) = pt';
            end
            eSeamXrot(iE,:,:) = R;
        end
    end
    
    % Part 2: Move extended edges
    for iE = 1:length(eSeamX)
        if eSeamX(iE)==1
            p1 = eSeamXptsRot(iE,[3:4]);
            p2lst = o(1).vOut(o(1).eOut(:,1),:);
            rho = sum((ones(size(p2lst,1),1)*p1 - p2lst).^2,2).^0.5;
            [~,idx] = min(rho);
            
            eSeamXptsRot(iE,[5:6]) = p2lst(idx,:);
            p0 = p2lst(idx,:) - p1;
            theta = eSeamXtheta(iE);
            R = [cos(theta) sin(theta); -sin(theta) cos(theta)];
            eSeamXpts(iE,[5:6]) = R * p0' + eSeamXpts(iE,3:4)';
        end
    end
    
    % Part 3: Assign fully populated matrices to the output struct
    lst = find(eSeamX==1);
    o(iO).eSeamExt = eSeamXpts(lst,[1 2 5 6]); 
    o(iO).eSeamExtRot = eSeamXptsRot(lst,[1 2 5 6]);
    
    o(iO).eSeamX = eSeamX;
    o(iO).eSeamXidx = eSeamXidx;
    o(iO).eSeamXpts = eSeamXpts;
    o(iO).eSeamXlen = eSeamXlen;
    o(iO).eSeamXtheta = eSeamXtheta;
    o(iO).eSeamXrot = eSeamXrot;
    o(iO).eSeamXptsRot = eSeamXptsRot;
end

%% --- Final Processing Steps ---
% Remove duplicates nodes and edges
for iO = 1:3
    if isempty(o(iO).v), continue; end
    v = o(iO).v;
    e = o(iO).e;
    [~, I, J] = unique(v, 'rows', 'first');
    map = J;
    edgesNew = map(e);
    edgesNew = edgesNew(edgesNew(:,1) ~= edgesNew(:,2), :);
    nodesNew = v(I,:);
    edgesNew = sort(edgesNew, 2);
    [~, I_e, ~] = unique(edgesNew, 'rows', 'first');
    o(iO).v = nodesNew;
    o(iO).e = edgesNew(I_e,:);
end

% Final Extension Calculations
for iO = 2:3
    if ~isfield(o(iO), 'eSeamExt') || isempty(o(iO).eSeamExt), continue; end
    v = o(iO).v;
    e = o(iO).e;
    for pt = 1:size(o(iO).eSeamExt,1)
        vidx = find(ismember(v, o(iO).eSeamExt(pt,1:2), 'rows'));
        if isempty(vidx), continue; end
        eidx = find(e(:,1) == vidx | e(:,2) == vidx);
        pts = [];
        vidxs = e(eidx,:);
        other_vidx = setdiff(vidxs(:),vidx);
        for u = 1:length(other_vidx)
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

if isfield(o(1), 'v')
    v = o(1).v;
    e = o(1).e;
    for iO = 2:3
        if ~isfield(o(iO), 'eSeamExtRot') || isempty(o(iO).eSeamExtRot), continue; end
        for pt = 1:size(o(iO).eSeamExtRot,1)
             if size(o(iO).eSeamExtRot,2) < 4, continue; end
            vidx = find(ismember(v, o(iO).eSeamExtRot(pt,3:4), 'rows'));
            if isempty(vidx), continue; end
            eidx = find(e(:,1) == vidx | e(:,2) == vidx);
            pts = [];
            vidxs = e(eidx,:);
            other_vidx = setdiff(vidxs(:),vidx);
            for u = 1:length(other_vidx)
                pts = [pts v(vidx,:) v(other_vidx(u),:)];
            end
            if ~isfield(o(iO), 'eExtensionsTop'), o(iO).eExtensionsTop(pt).pts = []; end
            o(iO).eExtensionsTop(pt).pts = pts;
            FixedPts = [o(iO).eSeamExt(pt,1) o(iO).eSeamExt(pt,2); o(iO).eSeamExt(pt,3) o(iO).eSeamExt(pt,4)];
            Movingpts = [o(iO).eSeamExtRot(pt,1) o(iO).eSeamExtRot(pt,2); o(iO).eSeamExtRot(pt,3) o(iO).eSeamExtRot(pt,4)];
            tform = fitgeotrans(Movingpts,FixedPts,'nonreflectivesimilarity');
            if ~isfield(o(iO).eExtensionsTop(pt), 'pts') || isempty(o(iO).eExtensionsTop(pt).pts), continue; end
            Xpts = o(iO).eExtensionsTop(pt).pts(1:2:end);
            Ypts = o(iO).eExtensionsTop(pt).pts(2:2:end);
            [Txpts,Typts] = transformPointsForward(tform, Xpts, Ypts);
            temppts = [Txpts;Typts];
            o(iO).eExtensionsTopRot(pt).pts = temppts(:)';
        end
    end
end


end
% --- END OF MAIN FUNCTION ---

% =========================================================================
% --- ALL HELPER FUNCTIONS ARE PLACED BELOW AS LOCAL FUNCTIONS ---
% =========================================================================

function [V, E] = connectionsToVE(connections)
    valid_coords = connections(~any(isnan(connections), 2), :);
    if isempty(valid_coords)
        V = [];
        E = [];
        return;
    end
    [V, ~, ic] = unique(valid_coords, 'rows', 'stable');
    E = reshape(ic, 2, [])';
end

function connections = find_connections(centers, positions, k, max_distance)
    connections = []; 
    if isempty(centers) || isempty(positions), return; end
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

% % CORRECTED a syntax error in this function definition
% function [detectors_left_side, sources_left_side, detectors_right_side, ...
%     sources_right_side, det_centroids_left, src_centroids_left, ...
%     det_centroids_right, src_centroids_right, ~, ~] = ...
%     get_UHD_optode_pos(sideLeftOutline, sideRightOutline, ~)
%     
%     vertical_spacing = 12;
%     outer_curve_factor = 1;
%     outer_curve_offset = outer_curve_factor*vertical_spacing;
%     for u = 1:2
%         if  u == 1
%             horizontal_spacing = vertical_spacing * sind(60);
%             curve_y = sideLeftOutline(:,2); 
%             curve_x = sideLeftOutline(:,1);
%             [outer_curve_x, outer_curve_y] = create_outer_curve(curve_x, curve_y, outer_curve_offset);
%         else
%             horizontal_spacing = -vertical_spacing * sind(60);
%             curve_y = sideRightOutline(:,2); 
%             curve_x = sideRightOutline(:,1);
%             [outer_curve_x, outer_curve_y] = create_outer_curve(curve_x, curve_y, -outer_curve_offset);
%         end
% 
%         detectors = [];
%         sources = [];
%         
%         polygon_x = [curve_x(1); curve_x; curve_x(end)];
%         polygon_y = [curve_y(1); curve_y; curve_y(end)];
% 
%         x_pos = curve_x(end); 
%         y_start_col1 = []; 
% 
%         for y_pos = max(curve_y)+outer_curve_factor*vertical_spacing:-vertical_spacing:min(curve_y)-outer_curve_factor*vertical_spacing
%             if isempty(y_start_col1)
%                y_start_col1 = y_pos;
%             end
%         end
%         
%         column_is_detector_only = false; 
%         x_pos = x_pos + horizontal_spacing;
% 
%         if u == 1
%             condition = @(x_pos) x_pos < max(outer_curve_x);
%         else
%             condition = @(x_pos) x_pos > min(outer_curve_x);
%         end
%         det_or_src = 0;
%         while condition(x_pos)
%             if column_is_detector_only
%                 for y_pos = max(curve_y)+outer_curve_factor*vertical_spacing:-vertical_spacing:min(curve_y)-outer_curve_factor*vertical_spacing
%                     if inpolygon(x_pos, y_pos, polygon_x, polygon_y)
%                         detectors = [detectors; x_pos, y_pos];
%                     end
%                 end
%             else
%                 y_start_col2 = y_start_col1 +(outer_curve_factor-0.5)*vertical_spacing;
%                 for i = 0:25
%                     if mod(outer_curve_factor+det_or_src,2) == 1
%                         detector_y = y_start_col2 - i * (2 * vertical_spacing);
%                         source_y = detector_y - vertical_spacing;
%                     else
%                         source_y = y_start_col2 - i * (2 * vertical_spacing);
%                         detector_y = source_y - vertical_spacing;
%                     end
%                     if inpolygon(x_pos, detector_y, polygon_x, polygon_y)
%                         detectors = [detectors; x_pos, detector_y]; 
%                     end
%                     if inpolygon(x_pos, source_y, polygon_x, polygon_y)
%                         sources = [sources; x_pos, source_y]; 
%                     end
%                 end
%                 det_or_src = det_or_src+1;
%             end
%             x_pos = x_pos + horizontal_spacing;
%             column_is_detector_only = ~column_is_detector_only;
%         end
%         
%         if u == 1
%             detectors_left_side = detectors;
%             sources_left_side = sources;
%         else
%             detectors_right_side = detectors;
%             sources_right_side = sources;
%         end
%     end
%     
%     % Placeholder centroid calculation
%     if ~isempty(sources_left_side)
%         det_centroids_left = mean(sources_left_side, 1) - [10.39, 6];
%         src_centroids_left = mean(sources_left_side, 1) + [10.39, -6];
%     else
%         det_centroids_left = [];
%         src_centroids_left = [];
%     end
%     if ~isempty(sources_right_side)
%         det_centroids_right = mean(sources_right_side, 1) + [10.39, 6];
%         src_centroids_right = mean(sources_right_side, 1) - [10.39, -6];
%     else
%         det_centroids_right = [];
%         src_centroids_right = [];
%     end
% end

function [outer_curve_x, outer_curve_y] = create_outer_curve(curve_x, curve_y, outer_curve_offset)
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