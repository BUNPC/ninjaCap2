function exclusionMask = createExclusionMask(v_conn, e_conn, maskSize, radius)
%createExclusionMask Creates a dilated mask, filling the interior of the structure.

initialMask = false(maskSize);

% Mark the vertices and edges on the mask (same as before)
if ~isempty(v_conn)
    v_round = round(v_conn);
    v_round(:,1) = max(min(v_round(:,1), maskSize(2)), 1);
    v_round(:,2) = max(min(v_round(:,2), maskSize(1)), 1);
    indices = sub2ind(maskSize, v_round(:,2), v_round(:,1));
    initialMask(indices) = true;
end
if ~isempty(e_conn)
    v_round = round(v_conn);
    for i = 1:size(e_conn, 1)
        p1 = v_round(e_conn(i, 1), :); p2 = v_round(e_conn(i, 2), :);
        x1 = p1(1); y1 = p1(2); x2 = p2(1); y2 = p2(2);
        dx = abs(x2 - x1); dy = -abs(y2 - y1);
        sx = sign(x2 - x1); sy = sign(y2 - y1); err = dx + dy;
        while true
            if (y1 >= 1 && y1 <= maskSize(1) && x1 >= 1 && x1 <= maskSize(2))
                initialMask(y1, x1) = true;
            end
            if (x1 == x2 && y1 == y2), break; end
            e2 = 2 * err;
            if (e2 >= dy), err = err + dy; x1 = x1 + sx; end
            if (e2 <= dx), err = err + dx; y1 = y1 + sy; end
        end
    end
end

% --- NEW: Fill the interior of the v_conn structure ---
if size(v_conn, 1) > 2
    try
        k = convhull(v_conn);
        interiorMask = poly2mask(v_conn(k,1), v_conn(k,2), maskSize(1), maskSize(2));
        initialMask = initialMask | interiorMask;
    catch
        disp('Could not generate convex hull; proceeding without filling interior.');
    end
end

% Dilate the mask to create the final exclusion zone
if radius > 0
    se = strel('disk', radius);
    exclusionMask = imdilate(initialMask, se);
else
    exclusionMask = initialMask;
end

end