function [transformed_main_curve, transformed_src_pts, transformed_det_pts] = curv_transform(main_curve_pts, target_line_pts, src_pts, det_pts)
    % The points in the data files may not be in order. We need to sort them
    % to trace the lines from end to end.

    disp('Ordering points and parameterizing curves...');

    % Order the points along each line
    main_curve_pts_ordered = order_points_along_line(main_curve_pts);
    target_line_pts_ordered = order_points_along_line(target_line_pts);

    % Calculate the cumulative distance (arc length) for each point
    main_curve_arclength = calculate_arc_lengths(main_curve_pts_ordered);
    target_line_arclength = calculate_arc_lengths(target_line_pts_ordered);

    % Scale the target line's parameterization to match the source curve's length.
    % This handles cases where the lines have slightly different total lengths.
    target_line_arclength = target_line_arclength * (main_curve_arclength(end) / target_line_arclength(end));


    %% --- Create the Transformation Map ---
    % We create an interpolant that can find the (x, y) coordinates on the
    % target line for any given arc length value from the source curve.

    disp('Creating transformation map...');

    % Create interpolation functions. Given an arc length 's', these will return
    % the corresponding X and Y coordinates on the target line.
    interp_target_x = @(s) interp1(target_line_arclength, target_line_pts_ordered(:,1), s, 'linear', 'extrap');
    interp_target_y = @(s) interp1(target_line_arclength, target_line_pts_ordered(:,2), s, 'linear', 'extrap');

    %% --- Transform the Main Curve and Other Rows ---
    disp('Applying transformation...');

    % a) Transform the main curve itself
    transformed_main_curve = zeros(size(main_curve_pts_ordered));
    for i = 1:size(main_curve_pts_ordered, 1)
        s = main_curve_arclength(i); % Get arc length of the current point
        transformed_main_curve(i, 1) = interp_target_x(s);
        transformed_main_curve(i, 2) = interp_target_y(s);
    end

    % b) Transform the adjacent rows
    % This is the most complex part. For each point in a side row, we find its
    % relationship to the main curve and re-apply it to the new straight line.
    transformed_src_pts = transform_adjacent_row(src_pts, main_curve_pts_ordered, transformed_main_curve);
    transformed_det_pts = transform_adjacent_row(det_pts, main_curve_pts_ordered, transformed_main_curve);


    %% --- Visualization ---
%     disp('Plotting results...');
% 
%     % Plot 1: The original configuration
%     figure('Name', 'Original Optode Configuration', 'Position', [100, 400, 800, 400]);
%     hold on;
%     plot(target_line_pts(:,1), target_line_pts(:,2), 'k-s', 'LineWidth', 2, 'MarkerFaceColor', 'k', 'DisplayName', 'Target Line');
%     plot(main_curve_pts(:,1), main_curve_pts(:,2), 'o-r', 'LineWidth', 1.5, 'MarkerFaceColor', 'r', 'DisplayName', 'Main Curve');
%     plot(row1_pts(:,1), row1_pts(:,2), 'd-b', 'LineWidth', 1, 'MarkerFaceColor', 'b', 'DisplayName', 'Row 1');
%     plot(row2_pts(:,1), row2_pts(:,2), '^-c', 'LineWidth', 1, 'MarkerFaceColor', 'c', 'DisplayName', 'Row 2');
%     title('Original Configuration');
%     xlabel('X Position'); ylabel('Y Position');
%     legend('show');
%     axis equal; grid on;
% 
%     % Plot 2: The transformed configuration
%     figure('Name', 'Transformed Optode Configuration', 'Position', [950, 400, 800, 400]);
%     hold on;
%     plot(target_line_pts(:,1), target_line_pts(:,2), 'k-s', 'LineWidth', 2, 'MarkerFaceColor', 'k', 'DisplayName', 'Target Line');
%     plot(transformed_main_curve(:,1), transformed_main_curve(:,2), 'o-r', 'LineWidth', 1.5, 'MarkerFaceColor', 'r', 'DisplayName', 'Transformed Main Curve');
%     plot(transformed_row1_pts(:,1), transformed_row1_pts(:,2), 'd-b', 'LineWidth', 1, 'MarkerFaceColor', 'b', 'DisplayName', 'Transformed Row 1');
%     plot(transformed_row2_pts(:,1), transformed_row2_pts(:,2), '^-c', 'LineWidth', 1, 'MarkerFaceColor', 'c', 'DisplayName', 'Transformed Row 2');
% 
%     % Optional: Draw lines showing the mapping from old to new positions
%     for k=1:size(row1_pts,1)
%         plot([row1_pts(k,1), transformed_row1_pts(k,1)], [row1_pts(k,2), transformed_row1_pts(k,2)], ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
%     end
% 
%     title('Transformed Configuration');
%     xlabel('X Position'); ylabel('Y Position');
%     legend('show');
%     axis equal; grid on;
end


%% --- Helper Functions ---

function ordered_pts = order_points_along_line(pts)
    % Sorts a cloud of 2D points to form a contiguous line by repeatedly
    % finding the nearest neighbor.
    num_pts = size(pts, 1);
    remaining_indices = 1:num_pts;
    ordered_indices = zeros(1, num_pts);
    
    % Start with the point that has the lowest Y, then X (a consistent start)
    [~, start_idx_in_all] = min(pts(:,2));
    ordered_indices(1) = start_idx_in_all;
    remaining_indices(start_idx_in_all) = [];
    
    for i = 2:num_pts
        last_pt_idx = ordered_indices(i-1);
        last_pt_coords = pts(last_pt_idx, :);
        
        % Find the nearest point among the remaining ones
        [~, nearest_idx_in_remaining] = min(sum((pts(remaining_indices, :) - last_pt_coords).^2, 2));
        
        % Get its original index
        next_pt_idx = remaining_indices(nearest_idx_in_remaining);
        
        ordered_indices(i) = next_pt_idx;
        remaining_indices(remaining_indices == next_pt_idx) = [];
    end
    
    ordered_pts = pts(ordered_indices, :);
end

function arc_lengths = calculate_arc_lengths(ordered_pts)
    % Calculates the cumulative distance along a path of ordered points.
    diffs = diff(ordered_pts, 1, 1);
    segment_lengths = vecnorm(diffs, 2, 2);
    arc_lengths = [0; cumsum(segment_lengths)];
end

function transformed_row = transform_adjacent_row(row_pts, source_curve, transformed_curve)
    % Transforms an adjacent row of points by preserving its offset from the main curve.
    transformed_row = zeros(size(row_pts));
    
    for i = 1:size(row_pts, 1)
        pt_Q = row_pts(i, :);
        
        % 1. Find the closest point P on the original source curve
        [~, closest_idx] = min(sum((source_curve - pt_Q).^2, 2));
        pt_P_source = source_curve(closest_idx, :);
        
        % 2. Calculate the original offset vector
        offset_vector = pt_Q - pt_P_source;
        
        % 3. Find the corresponding point on the transformed curve
        pt_P_transformed = transformed_curve(closest_idx, :);
        
        % 4. Calculate the rotation needed for the offset vector
        % Get the local direction (tangent) of the original curve at P
        tangent_source = get_tangent(source_curve, closest_idx);
        
        % Get the local direction (tangent) of the new curve at P'
        tangent_transformed = get_tangent(transformed_curve, closest_idx);
        
        % Find the angle of each tangent vector
        angle_source = atan2(tangent_source(2), tangent_source(1));
        angle_transformed = atan2(tangent_transformed(2), tangent_transformed(1));
        
        % The angle by which we need to rotate the offset
        rotation_angle = angle_transformed - angle_source;
        
        % 5. Rotate the offset vector
        R = [cos(rotation_angle), -sin(rotation_angle); 
             sin(rotation_angle),  cos(rotation_angle)];
        rotated_offset = (R * offset_vector')';
        
        % 6. The new point is the transformed curve point plus the rotated offset
        transformed_row(i, :) = pt_P_transformed + rotated_offset;
    end
end

function tangent = get_tangent(curve_pts, idx)
    % Calculates the local tangent at a point on a curve using finite differences.
    num_pts = size(curve_pts, 1);
    if idx == 1 % Forward difference for the start point
        tangent = curve_pts(2, :) - curve_pts(1, :);
    elseif idx == num_pts % Backward difference for the end point
        tangent = curve_pts(end, :) - curve_pts(end-1, :);
    else % Central difference for interior points
        tangent = curve_pts(idx+1, :) - curve_pts(idx-1, :);
    end
    tangent = tangent / norm(tangent); % Normalize to a unit vector
end
