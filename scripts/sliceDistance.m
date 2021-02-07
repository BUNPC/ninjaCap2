function distance = sliceDistance(slice, p1, p2)
%SLICEDISTANCE Summary of this function goes here
%   Detailed explanation goes here

% convert to connectivity
[slice_vertices, ~, slice_connect] = edgesToConnectivity(slice);

% opening vertex
[~, v_start] = min(sum(bsxfun(@minus, p1, slice_vertices) .^ 2, 2));
[~, v_end] = min(sum(bsxfun(@minus, p2, slice_vertices) .^ 2, 2));

% sets
set_closed = false(1, size(slice_vertices, 1));
set_open = false(1, size(slice_vertices, 1));
set_open(v_start) = true;

% scores
score_g = inf(1, size(slice_vertices, 1));
score_g(v_start) = 0;
score_f = inf(1, size(slice_vertices, 1));
score_f(v_start) = heuristic(slice_vertices(v_start, :), slice_vertices(v_end, :));

% allow reconstructing path (not required for distance, but good for
% debugging)
came_from = zeros(1, size(slice_vertices, 1));

while any(set_open)
    % find lowest f score
    o = find(set_open);
    [~, idx] = min(score_f(set_open));
    current = o(idx);
    if current == v_end
        break;
    end
    
    % remove from open
    set_open(current) = false;
    set_closed(current) = true;
    
    % for neighbors
    for neighbor = find(slice_connect(current, :))
        % ignore already visited points
        if set_closed(neighbor)
            continue
        end
        
        % calculate g score (distance to current + distance from current to
        % neighbor)
        g = score_g(current) + sqrt(sum((slice_vertices(current, :) - slice_vertices(neighbor, :)) .^ 2, 2));
        
        if ~set_open(neighbor)
            set_open(neighbor) = true;
        elseif g >= score_g(neighbor)
            % no improvement
            continue;
        end
        
        % update f and g scores
        score_g(neighbor) = g;
        score_f(neighbor) = g + heuristic(slice_vertices(neighbor, :), slice_vertices(v_end, :));
        came_from(neighbor) = current;
    end
end

% debug: reconstruct path
% current = v_end;
% path = [];
% figure;
% plot(slice(:, 1), slice(:, 2));
% hold on;
% while came_from(current) ~= 0
%     path = [path; came_from(current) current];
%     plot(slice_vertices([came_from(current); current], 1), slice_vertices([came_from(current); current], 2), 'r');
%     current = came_from(current);
% end
% hold off;

% distance from v_start to v_end
distance = score_g(v_end);

    function h = heuristic(pt, dest)
        h = sqrt(sum((pt - dest) .^ 2, 2));
    end
end

