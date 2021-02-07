function vertices = eliminateShortEdges(vertices, idx, threshold)
%ELIMINATESHORTEDGES Summary of this function goes here
%   Detailed explanation goes here

num_vertices = size(vertices, 1);

% make sure it is a lattice of line segments
assert(0 == mod(num_vertices, 3));
assert(all(all(isnan(vertices(3:3:end, :)))));

% make sure idx matches vertices
assert(length(idx) == num_vertices);

% find which ones are edges
is_edge = idx(1:3:end) | idx(2:3:end);

% calculate lengths (waste some energy calculating all lengths, not just
% edge pieces, but small beans)
lengths = sqrt(sum((vertices(1:3:end, :) - vertices(2:3:end, :)) .^ 2, 2));

% short edges
to_eliminate = is_edge & lengths <= threshold;

for i = find(to_eliminate)'
    % get vertices
    v1 = (i - 1) * 3 + 1;
    v2 = (i - 1) * 3 + 2;
    
    % make v2 be the edge point (the one to use when updating all other
    % edges)
    
    if idx(v1) && idx(v2)
        % if both are edges (super short line, consider the one that has
        % less vertices to be the edge point)
        cnt1 = sum(ismember(vertices, vertices(v1, :), 'rows'));
        cnt2 = sum(ismember(vertices, vertices(v2, :), 'rows'));
        
        % swap them so v2 has fewer vertices
        if cnt2 > cnt1
            [v2, v1] = deal(v1, v2);
        end
    elseif idx(v1)
        % swap them so v2 has the edge point
        [v2, v1] = deal(v1, v2);
    end
    
    % find entries to replace
    to_replace = sum(bsxfun(@minus, vertices, vertices(v1, :)) .^ 2, 2) < eps;
    vertices(to_replace, :) = repmat(vertices(v2, :), sum(to_replace), 1);
end

% eliminate short lines
vertices(repelem(to_eliminate, 3), :) = [];

end
