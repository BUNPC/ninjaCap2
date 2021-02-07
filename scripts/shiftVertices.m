function vertices = shiftVertices(vertices, dist, desired_dist, line, reinforce_threshold, fully_connected)
%SHIFTVERTICES Moves vertices for aligned stitching
%   Vertices are all vertices in a grid. For those along the line in
%   question (`line`), vector `dist` contains the position along the line
%   (and should be nan for vertices not along line). Vertices along the
%   line are shifted so that the distances align with `desired_dist`.

%% validate inputs

num_vertices = size(vertices, 1);

% make sure it is a lattice of line segments
assert(0 == mod(num_vertices, 3));
assert(all(all(isnan(vertices(3:3:end, :)))));

% make sure idx matches vertices
assert(length(dist) == num_vertices);

% make sure the line segment is defined
assert(2 == size(line, 2));

% default
if ~exist('reinforce_threshold', 'var')
    reinforce_threshold = [];
end

%% prep work

% calculate distances for each line segment (used to map distances back to
% line)
lines_num = size(line, 1) - 1;
lines_dist = zeros(lines_num, 1);
for i = 1:lines_num
    lines_dist(i) = sqrt(sum((line(i + 1, :) - line(i, :)) .^ 2));
end

% generate new vertices
vertices_new = zeros(length(desired_dist), 2);

% sort desired distances
desired_dist = sort(desired_dist);

d = 0;
l = 0;
for i = 1:length(desired_dist)
    % move along line to find segment that will contain desired distance
    % (depends on the fact that desired_dist are sorted, so we know that
    % the current desired distance must be on or past the previous one)
    while desired_dist(i) > d && l < lines_num
        l = l + 1;
        d = d + lines_dist(l);
    end
    
    % figure out fraction
    frac = max(d - desired_dist(i), 0) / lines_dist(l);
    
    % calculate vertex
    vertices_new(i, :) = frac * line(l, :) + (1 - frac) * line(l + 1, :);
end

%% actually map distances to nearest points

% calculate distance matrix
dm = abs(dist(:) - desired_dist(:)');

% find closest vertex
[c_d, closest] = min(dm, [], 2);

% vertices to move
idx = ~isnan(dist);

% move vertices
vertices(idx, :) = vertices_new(closest(idx), :);

%% fully connected

if fully_connected
    connections = accumarray(closest(idx), 1, [size(desired_dist, 1), 1]);
    unconnected = find(connections == 0);
    for i = 1:length(unconnected)
        idx = unconnected(i);
        
        % find nearby connected points
        nearby = [find(connections(1:idx), 1, 'last'); ...
            idx - 1 + find(connections(idx:end), 1)];
        
        % nothing can be done
        if isempty(nearby)
            continue;
        end
        
        % nearby vertices (based on the previous closest point data)
        nearby_idx = find(ismember(closest, nearby));
        
        % switch to the other end of the lines (if pointing to start of
        % line [x % 3 == 1], then add 1; if pointing to end of line [x % 3
        % == 2], then subtract one)
        s = mod(nearby_idx, 3); s(s == 2) = -1;
        nearby_idx = nearby_idx + s;
        
        % find which one is closest
        nearby_dist = sum(bsxfun(@minus, vertices(nearby_idx, :), vertices_new(idx, :)) .^ 2, 2);
        [~, to_connect] = min(nearby_dist);
        to_connect = nearby_idx(to_connect);
        
        % add connection
        vertices = [vertices; vertices(to_connect, :); vertices_new(idx, :); nan nan]; %#ok<AGROW>
    end
end

%% add another support

% if a big shift, add vertices to other nearest point
if ~isempty(reinforce_threshold)
    for i = find(c_d > reinforce_threshold & idx)'
        if mod(i, 3) == 2
            v1 = vertices(i - 1, :);
        else
            v1 = vertices(i + 1, :);
        end

        if desired_dist(closest(i)) < dist(i)
            if closest(i) == length(desired_dist)
                continue;
            end
            v2 = vertices_new(closest(i) + 1, :);
        else
            if closest(i) == 1
                continue;
            end
            v2 = vertices_new(closest(i) - 1, :);
        end

        % add vertex
        vertices = [vertices; v1; v2; nan nan]; %#ok<AGROW>
    end
end

end

