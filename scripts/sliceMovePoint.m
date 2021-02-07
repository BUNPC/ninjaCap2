function p1 = sliceMovePoint(slice, p0, rel_to, dist)
%SLICEMOVEPOINT Summary of this function goes here
%   Detailed explanation goes here

% convert to connectivity
[vertices, ~, connectivity] = edgesToConnectivity(slice);

% find closest point
[v0, i0] = closestVertex(vertices, p0);

% visited
visited = false(1, size(vertices, 1));

% find options
potentials = expandVertex(i0, 0);
if isempty(potentials)
    error('Unable to move point.');
end

% get distances
distances = sqrt(sum(bsxfun(@minus, potentials, rel_to) .^ 2, 2));
distances_rel = distances - norm(v0 - rel_to, 2);

% sort 'em
[distances_rel, idx] = sort(distances_rel);

% allow moving towards or away
if sign(dist) == 1 % positive means move further away
    if distances_rel(end) < 0
        error('Only found points closer.');
    end
    
    p1 = potentials(idx(end), :);
else
    if distances_rel(1) > 0
        error('Only found points further.');
    end
    
    p1 = potentials(idx(1), :);
end

    function res = expandVertex(i, d)
        % mark as visited
        visited(i) = true;
        
        % trvaeled far enough
        if d > abs(dist)
            res = vertices(i, :);
            return;
        end
        
        % find connections (that have not been visited
        res = [];
        other_is = find(connectivity(i, :) & ~visited);
        for other_i = other_is
            new_d = d + norm(vertices(i, :) - vertices(other_i, :), 2);
            res = [res; expandVertex(other_i, new_d)];
        end
    end

end

