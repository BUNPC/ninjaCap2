function [detectors_left_side, sources_left_side, detectors_right_side, ...
    sources_right_side, det_centroids_left, src_centroids_left, ...
    det_centroids_right, src_centroids_right, outer_sources_left_side, ...
    outer_sources_right_side, outer_detectors_left_side, outer_detectors_right_side] = ...
    get_UHD_optode_pos_modified(sideLeftOutline, sideRightOutline)
%
% This is the fully corrected version of the function.
%
    %% ---Fill optodes on side Left ---
    vertical_spacing = 12;
    outer_curve_factor = 1;
    outer_curve_offset = outer_curve_factor*vertical_spacing;

    % Initialize output variables
    outer_detectors_left_side = []; outer_sources_left_side = [];
    outer_detectors_right_side = []; outer_sources_right_side = [];

    for u = 1:2
        if  u == 1
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

        % ** CORRECTED ASSIGNMENT BLOCK **
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
    end % End of the for loop

    %% Find centroids of detector triangle groups
    det_centroids_left = get_det_centroids(detectors_left_side, 12, 15);
    det_centroids_right = get_det_centroids(detectors_right_side, 12, 15);

    src_centroids_left = det_centroids_left + [12*sind(60) 0];
    src_centroids_right = det_centroids_right - [12*sind(60) 0];
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