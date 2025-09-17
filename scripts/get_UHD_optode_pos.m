function [detectors_left_side, sources_left_side, detectors_right_side, ...
    sources_right_side, det_centroids_left, src_centroids_left, ...
    det_centroids_right, src_centroids_right, outer_sources_left_side, outer_sources_right_side] = ...
    get_UHD_optode_pos(sideLeftOutline, sideRightOutline, topOutline)
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


        % Store coordinates for plotting
        detectors = [];
        sources = [];

        detectors_all = [];
        sources_all = [];


        %%
        % --- Setup the Plot and Boundary ---
        % Create a closed polygon for boundary checking using the 'inpolygon' function.
        polygon_x = [curve_x(1); curve_x; curve_x(end)];
        polygon_y = [curve_y(1); curve_y; curve_y(end)];

        outer_polygon_x = [outer_curve_x(1); outer_curve_x; outer_curve_x(end)];
        outer_polygon_y = [outer_curve_y(1); outer_curve_y; outer_curve_y(end)];
        %%

        % Start the first column exactly at the flat edge (x=0)
        x_pos = curve_x(end); 
        y_start_col1 = []; % To store the y-position of the first detector

        % Start from the top-most point of the curve and place detectors downwards
        for y_pos = max(curve_y)+outer_curve_factor*vertical_spacing:-vertical_spacing:min(curve_y)-outer_curve_factor*vertical_spacing
            % At x=0, all points are within the vertical bounds of the shape
            detectors = [detectors; x_pos, y_pos];
            if isempty(y_start_col1)
               y_start_col1 = y_pos; % Save the y-position of the topmost detector
            end
    %         if inpolygon(x_pos, y_pos, outer_polygon_x, outer_polygon_y)
               detectors_all = [detectors_all; x_pos, y_pos];
    %         end
        end
        detectors = detectors(outer_curve_factor+1:end-outer_curve_factor,:);

        % --- Place Subsequent Columns ---

        % Initialize for the loop. The first column was detectors, so the next is source/detector.
        column_is_detector_only = false; 
        x_pos = x_pos + horizontal_spacing;

        if u == 1
            condition = @(x_pos) x_pos < max(outer_curve_x);
        else
            condition = @(x_pos) x_pos > min(outer_curve_x);
        end
        det_or_src = 0;
        while condition(x_pos)
            % Continue until we pass the right-most point of the curve

            if column_is_detector_only
                % --- Place a column of DETECTORS ---

                % Loop through the entire y-range of the shape and let inpolygon check
                for y_pos = max(curve_y)+outer_curve_factor*vertical_spacing:-vertical_spacing:min(curve_y)-outer_curve_factor*vertical_spacing
                     if inpolygon(x_pos, y_pos, polygon_x, polygon_y)
                        detectors = [detectors; x_pos, y_pos];
                     end
                     if inpolygon(x_pos, y_pos, outer_polygon_x, outer_polygon_y) || inpolygon(x_pos, y_pos, polygon_x, polygon_y)
                       detectors_all = [detectors_all; x_pos, y_pos];
                    end
                end

            else
                % --- Place a column of SOURCE and DETECTOR pairs ---

                % The first detector's position is centered between the first two 
                % detectors of the first column.
                y_start_col2 = y_start_col1 +(outer_curve_factor-0.5)*vertical_spacing;

                % Place pairs starting from this y_start_col2 and moving outwards.

                % Loop to place pairs downwards from the starting point (including start)
                for i = 0:25 % Arbitrary limit, boundary check will stop it
                    if mod(outer_curve_factor+det_or_src,2) == 1
                        detector_y = y_start_col2 - i * (2 * vertical_spacing);
                        source_y = detector_y - vertical_spacing; % Source is 12mm below its detector
                    else
                        source_y = y_start_col2 - i * (2 * vertical_spacing);
                        detector_y = source_y - vertical_spacing; % Source is 12mm below its detector
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

            % --- Update for the next iteration ---
            x_pos = x_pos + horizontal_spacing; % Move to the next column
            column_is_detector_only = ~column_is_detector_only; % Alternate column type
        end
        nearby_src_pts = pts_closer_to_curve([curve_x curve_y], sources_all, 35);
        nearby_det_pts = pts_closer_to_curve([curve_x curve_y], detectors_all, 35);
        if  u == 1
            detectors_left_side = detectors;
            sources_left_side = sources;
            [sideLeftOutline_top, sources_left_side_top, detectors_left_side_top] = curv_transform([curve_x curve_y], topOutline(1:20,:), nearby_src_pts, nearby_det_pts);
            outer_detectors_left_side = detectors_all;
            outer_sources_left_side = sources_all;
            sideLeftOutline_outer = [outer_polygon_x outer_polygon_y];
        else
            detectors_right_side = detectors;
            sources_right_side = sources;
            [sideRightOutline_top, sources_right_side_top, detectors_right_side_top] = curv_transform([curve_x curve_y], topOutline(26:46,:), nearby_src_pts, nearby_det_pts);
            outer_detectors_right_side = detectors_all;
            outer_sources_right_side = sources_all;
            sideRightOutline_outer = [outer_polygon_x outer_polygon_y];
        end
    end
    %% Find centroids of detector triangle groups

    det_centroids_left = get_det_centroids(outer_detectors_left_side,12, 15);
    det_centroids_right = get_det_centroids(outer_detectors_right_side,12, 15);

    src_centroids_left = det_centroids_left+[13.833 0];  %FIXME update hardcoded 13.883 to formula
    src_centroids_right = det_centroids_right-[13.833 0];  %FIXME update hardcoded 13.883 to formula

    %%
%     figure;
%     hold on; 
% 
%     % Plot the boundaries
%     plot(sideLeftOutline(:,1), sideLeftOutline(:,2), 'k-', 'LineWidth', 1.5);
%     plot(sideRightOutline(:,1), sideRightOutline(:,2), 'k-', 'LineWidth', 1.5);
%     plot(topOutline(:,1), topOutline(:,2), 'k-', 'LineWidth', 1.5);
%     plot(sideLeftOutline_outer(:,1), sideLeftOutline_outer(:,2), 'k--', 'LineWidth', 0.5);
%     plot(sideRightOutline_outer(:,1), sideRightOutline_outer(:,2), 'k--', 'LineWidth', 0.5);
%     % Set plot properties
%     title('Detector and Source Placement');
%     axis equal; % Ensure the aspect ratio is 1:1
%     % if ~isempty(detectors_left_side)
%     %     plot(detectors_left_side(:,1), detectors_left_side(:,2), 'bo', 'MarkerFaceColor', 'b', 'DisplayName', 'Detectors');
%     % end
% 
%     if ~isempty(sources_left_side)
%         plot(sources_left_side(:,1), sources_left_side(:,2), 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Sources');
%     end
% 
%     % if ~isempty(detectors_right_side)
%     %     plot(detectors_right_side(:,1), detectors_right_side(:,2), 'bo', 'MarkerFaceColor', 'b', 'DisplayName', 'Detectors');
%     % end
% 
%     if ~isempty(sources)
%         plot(sources_right_side(:,1), sources_right_side(:,2), 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Sources');
%     end
% 
%     % if ~isempty(detectors_left_side)
%     %     plot(outer_detectors_left_side(:,1), outer_detectors_left_side(:,2), 'bo', 'MarkerFaceColor', 'b', 'DisplayName', 'Detectors');
%     % end
% 
%     if ~isempty(sources_left_side)
%         plot(outer_sources_left_side(:,1), outer_sources_left_side(:,2), 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Sources');
%     end
% 
%     % if ~isempty(detectors_right_side)
%     %     plot(outer_detectors_right_side(:,1), outer_detectors_right_side(:,2), 'bo', 'MarkerFaceColor', 'b', 'DisplayName', 'Detectors');
%     % end
% 
%     if ~isempty(sources)
%         plot(outer_sources_right_side(:,1), outer_sources_right_side(:,2), 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Sources');
%     end
% 
%     if ~isempty(det_centroids_left)
%         plot(det_centroids_left(:,1), det_centroids_left(:,2), 'b*', 'MarkerFaceColor', 'r', 'DisplayName', 'Centroids');
%     end
% 
%     if ~isempty(det_centroids_right)
%         plot(det_centroids_right(:,1), det_centroids_right(:,2), 'b*', 'MarkerFaceColor', 'r', 'DisplayName', 'Centroids');
%     end
% 
%     if ~isempty(src_centroids_left)
%         plot(src_centroids_left(:,1), src_centroids_left(:,2), 'r*', 'MarkerFaceColor', 'r', 'DisplayName', 'Centroids');
%     end
% 
% 
%     if ~isempty(src_centroids_right)
%         plot(src_centroids_right(:,1), src_centroids_right(:,2), 'r*', 'MarkerFaceColor', 'r', 'DisplayName', 'Centroids');
%     end
% 
% 
%     plot(sideLeftOutline_top(:,1), sideLeftOutline_top(:,2), 'g-', 'LineWidth', 1.5);
%     plot(sideRightOutline_top(:,1), sideRightOutline_top(:,2), 'g-', 'LineWidth', 1.5);
% 
%     % if ~isempty(detectors_left_side)
%     %     plot(detectors_left_side_top(:,1), detectors_left_side_top(:,2), 'b*', 'MarkerFaceColor', 'b', 'DisplayName', 'Detectors');
%     % end
%     % 
%     % if ~isempty(sources_left_side)
%     %     plot(sources_left_side_top(:,1), sources_left_side_top(:,2), 'r*', 'MarkerFaceColor', 'r', 'DisplayName', 'Sources');
%     % end
%     % 
%     % if ~isempty(detectors_right_side)
%     %     plot(detectors_right_side_top(:,1), detectors_right_side_top(:,2), 'b*', 'MarkerFaceColor', 'b', 'DisplayName', 'Detectors');
%     % end
%     % 
%     % if ~isempty(sources)
%     %     plot(sources_right_side_top(:,1), sources_right_side_top(:,2), 'r*', 'MarkerFaceColor', 'r', 'DisplayName', 'Sources');
%     % end
% 
%     hold off;
end

function nearby_pts = pts_closer_to_curve(curve_pts, pts, max_dist)
    num_candidate_pts = size(pts, 1);
    min_distances = zeros(num_candidate_pts, 1);
    num_segments = size(curve_pts, 1) - 1;
    for i = 1:num_candidate_pts
        current_pt = pts(i, :);

        % Calculate the distance to every segment on the main curve

        dist_to_all_segments = zeros(num_segments, 1);

        for j = 1:num_segments
            v1 = curve_pts(j, :);
            v2 = curve_pts(j+1, :);
            dist_to_all_segments(j) = point_to_segment_distance(current_pt, v1, v2);
        end

        % The distance to the curve is the minimum of these segment distances
        min_distances(i) = min(dist_to_all_segments);
    end
    is_nearby = min_distances <= max_dist;
    nearby_pts = pts(is_nearby, :);
end

function dist = point_to_segment_distance(pt, v1, v2)
    % Calculates the shortest distance from a point (pt) to a line segment (v1, v2).
    v1v2 = v2 - v1;
    ptv1 = pt - v1;

    L2 = sum(v1v2.^2); % Squared length of the segment
    if L2 == 0 % The segment is actually a point
        dist = norm(ptv1);
        return;
    end

    % Project the point onto the infinite line defined by the segment.
    % The parameter 't' indicates where the projection lies.
    % t=0 -> projection is at v1
    % t=1 -> projection is at v2
    % 0<t<1 -> projection is between v1 and v2
    t = dot(ptv1, v1v2) / L2;

    % Clamp t to the range [0, 1] to find the closest point on the segment
    t = max(0, min(1, t));

    % Calculate the coordinates of the closest point on the segment
    projection = v1 + t * v1v2;

    % The distance is the norm of the vector from the point to its projection
    dist = norm(pt - projection);
end

function [outer_curve_x, outer_curve_y] = create_outer_curve(curve_x, curve_y, outer_curve_offset)

    % Calculate the tangent vectors for each segment of the inner curve.
    % We use central differences for more accuracy on interior points.
    dx = gradient(curve_x);
    dy = gradient(curve_y);

    % Calculate the normal vectors by rotating the tangent vectors 90 degrees.
    % For a curve on the right, the outward-pointing normal is (dy, -dx).
    normal_x = dy;
    normal_y = -dx;

    % Normalize the normal vectors to get unit vectors
    magnitude = sqrt(normal_x.^2 + normal_y.^2);
    % To prevent division by zero if any points are duplicates
    magnitude(magnitude == 0) = 1; 
    unit_normal_x = normal_x ./ magnitude;
    unit_normal_y = normal_y ./ magnitude;

    % Create the outer curve by moving each point along its normal vector
    outer_curve_x = curve_x + unit_normal_x * outer_curve_offset;
    outer_curve_y = curve_y + unit_normal_y * outer_curve_offset;
end

